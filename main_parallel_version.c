/*
  This file is part of mpiBWA

  NGS aligner inspired by BWA

  The project was developped by Frederic Jarlier from Institut Curie and Nicolas Joly from Institut Pasteur

  Copyright (C) 2016-2017  Institut Curie / Institut Pasteur



  This program is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the Free
  Software Foundation, either version 3 of the License, or (at your option)
  any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#define _GNU_SOURCE

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <inttypes.h> /* For PRIu64 */
#include <libgen.h>
#include <limits.h>   /* For PATH_MAX */
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include "bwa.h"
#include "bwamem.h"
#include "utils.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif


/* We require 64bit offsets */
#ifndef MPI_OFFSET
#define MPI_OFFSET MPI_LONG_LONG
#endif

#ifdef TIMING
#define xfprintf fprintf
#else
#define xfprintf(...) /**/
#endif

/* Default read size */
#define READSIZE (10UL * 1024 * 1024)
#define NUM_PROC_OFFSET 2
#define MAGIC_UNLOCK_TAG 0xBEEF

#define MPI_ASSERT(fx) do { \
  assert(fx == MPI_SUCCESS); \
} while(0)

typedef struct rma_mutex
{
  int owner, rank, size;
  char *req, *req_slice_buffer;
  MPI_Win win;
  MPI_Comm comm;
  MPI_Datatype req_slice_type;
} rma_mutex_t;


/**
 * Non-public: Create the MPI derived type corresponding to the req array on
 * the owner rank, minus the element mapped to the origin rank (this is written
 * to the `req_slice_type` member of partially constructed mutex `new`).
 */
static void create_req_buffer_type_offset(rma_mutex_t *new)
{
  int ret, nblock, *len, *dsp;

  new->req_slice_type = MPI_DATATYPE_NULL;

  len = NULL;
  dsp = NULL;

  if (new->rank == 0 || new->rank == new->size - 1)
    nblock = 1;
  else
    nblock = 2;

  len = calloc(nblock, sizeof(int));
  assert(len != NULL);
  dsp = calloc(nblock, sizeof(int));
  assert(dsp != NULL);

  if (new->rank == 0) {
    len[0] = new->size - 1;
    dsp[0] = 1;
  } else if (new->rank == new->size - 1) {
    len[0] = new->size - 1;
    dsp[0] = 0;
  } else {
    len[0] = new->rank;
    dsp[0] = 0;
    len[1] = new->size - 1 - new->rank;
    dsp[1] = new->rank + 1;
  }

  MPI_ASSERT(MPI_Type_indexed(nblock, len, dsp, MPI_BYTE, &new->req_slice_type));
  MPI_ASSERT(MPI_Type_commit(&new->req_slice_type));

  free(len);
  free(dsp);
}

/**
 * Initialize the RMA mutex `m` owned by rank `owner` of communcator `comm`.
 * This call is collective.
 */
void rma_mutex_init_offset(rma_mutex_t *m, MPI_Comm comm, int owner)
{
  int ret;
  MPI_Aint req_size;
  rma_mutex_t new;

  new.owner = owner;
  new.comm = comm;

  new.req = NULL;
  new.req_slice_buffer = NULL;

  MPI_ASSERT(MPI_Comm_rank(new.comm, &new.rank));
  MPI_ASSERT(MPI_Comm_size(new.comm, &new.size));

  new.req_slice_buffer = malloc(new.size);
  assert(new.req_slice_buffer != NULL);

  create_req_buffer_type_offset(&new);

  if (new.rank == new.owner) {
    req_size = new.size;
    MPI_ASSERT(MPI_Alloc_mem(req_size, MPI_INFO_NULL, &new.req));
    bzero(new.req, req_size);
  } else {
    req_size = 0;
  }

  MPI_ASSERT(MPI_Win_create(new.req, req_size, 1, MPI_INFO_NULL, new.comm, &new.win));

  /* lookin good - safe to modify argument */
  memcpy(m, &new, sizeof(rma_mutex_t));
}

/**
 * Free resources associated with the RMA mutex `m`. This call is collective.
 */
void rma_mutex_free_offset(rma_mutex_t *m)
{
  int ret;

  MPI_ASSERT(MPI_Win_free(&m->win)); /* collective */

  if (m->rank == m->owner)
    MPI_ASSERT(MPI_Free_mem(m->req));

  MPI_ASSERT(MPI_Type_free(&m->req_slice_type));

  free(m->req_slice_buffer);
}

/**
 * Acquire the lock associated with RMA mutex `m`.
 */
void rma_mutex_lock_offset(rma_mutex_t *m)
{
  int i, contested, ret;
  uint8_t one = 1;

  /* exclusive lock
   * .. get full req_buffer
   * .. set my entry to 1 */
  MPI_ASSERT(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, m->owner, 0, m->win));
  MPI_ASSERT(MPI_Get(m->req_slice_buffer, m->size - 1, MPI_BYTE, m->owner, 0, 1, m->req_slice_type, m->win));
  MPI_ASSERT(MPI_Put(&one, 1, MPI_BYTE, m->owner, m->rank, 1, MPI_BYTE, m->win));
  MPI_ASSERT(MPI_Win_unlock(m->owner, m->win));

  /* lock obtained? */
  contested = 0;
  for (i = 0; i < m->size - 1; i++)
    if (m->req_slice_buffer[i]) {
      contested = 1;
      break;
    }

  /* no, it's contested - wait for it ... */
  if (contested)
    MPI_ASSERT(MPI_Recv(NULL, 0, MPI_BYTE, MPI_ANY_SOURCE, MAGIC_UNLOCK_TAG, m->comm, MPI_STATUS_IGNORE));
}

/**
 * Non-public: Once the `req_slice_buffer` field of mutex `m` has been
 * populated, determine the next rank to receive the lock.
 * Return value: rank on the next process to get the lock (m->rank if no rank is
 * currently waiting).
 */
static int next_rank_offset(rma_mutex_t *m)
{
  int i, j, next;
  next = m->rank;
  for (i = 0; i < m->size - 1; i++) {
    j = (m->rank + i) % (m->size - 1);
    if (m->req_slice_buffer[j]) {
      if (j < m->rank)
        next = j;
      else
        next = j + 1;
      break;
    }
  }
  return next;
}

/**
 * Release the lock associated with RMA mutex `m`.
 */
void rma_mutex_unlock_offset(rma_mutex_t *m)
{
  int next, ret;
  uint8_t zero = 0;

  /* exclusive lock
   * .. get new full req_buffer
   * .. set my entry to 0 */
  MPI_ASSERT(MPI_Win_lock(MPI_LOCK_EXCLUSIVE, m->owner, 0, m->win));
  MPI_ASSERT(MPI_Get(m->req_slice_buffer, m->size - 1, MPI_BYTE, m->owner, 0, 1, m->req_slice_type, m->win));
  MPI_ASSERT(MPI_Put(&zero, 1, MPI_BYTE, m->owner, m->rank, 1, MPI_BYTE, m->win));
  MPI_ASSERT(MPI_Win_unlock(m->owner, m->win));

  /* who gets it next? */
  next = next_rank_offset(m);

  /* someone ... give it to them ... */
  if (next != m->rank)
    MPI_ASSERT(MPI_Send(NULL, 0, MPI_BYTE, next, MAGIC_UNLOCK_TAG, m->comm));
}


int main(int argc, char *argv[]) {
	double bef, aft;

	char *file_r1 = NULL, *file_r2 = NULL;
	char *buffer_r1, *buffer_r2;
	buffer_r1 = buffer_r2 = NULL;
	char *buffer_out     = NULL;
	MPI_File fd_tmp, fd_tmp1;
	uint8_t *a, *addr, *addr_map;
	size_t rlen_r1, rlen_r2;

	char *file_out       = NULL;
	size_t buffer_out_sz = 0;
	char *file_ref = NULL;

	char file_map[PATH_MAX], file_tmp[PATH_MAX];
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	off_t val_succ1;
	off_t val_succ2;
	int successor = 0;
	int name_len;
	off_t val1 = 0;
	off_t val2 = 0;
	int zero = 0;
	size_t curr[6];
	size_t csiz = sizeof(curr);
	ssize_t ssiz;
	size_t size_r1, size_r2, line_r1, line_r2, line_number;
	int count_r1, count_r2, full_r1, full_r2;

	int proc_num, rank_num, rank_shr;
	size_t localsize;
	size_t n = 0;
	struct flock lck;

	MPI_Aint size_shr;
	MPI_Comm comm_shr;

	MPI_File fh_map;
	MPI_File fh_r1, fh_r2, fh_r3, fh_r4, fh_out;

	MPI_Offset m, size_map, size_tot;
	MPI_Status status;
	MPI_Win win_shr = NULL;
	mem_opt_t *opt, opt0;
	mem_pestat_t pes[4], *pes0 = NULL;
	bwaidx_t indix;
	bseq1_t *seqs;

	int res, count;
	int files, nargs;
	size_t total_writen_by_proc = 0;
	size_t total_writen = 0;
	char *p, *q, *e;

	int c, copy_comment = 0;
	char *rg_line = NULL, *hdr_line = NULL;
	int ignore_alt = 0;
	const char *mode = NULL;

	rma_mutex_t mutex_offset;

	char *progname = basename(argv[0]);

	if (argc < 2) {
		fprintf(stderr, "Program: MPI version of BWA MEM\n\n"
			"Version: v%s\n\n"
			"Contact 1: Frederic Jarlier (frederic.jarlier@curie.fr)\n\n"
			"Contact 2: Nicolas Joly (njoly@pasteur.fr)\n\n"
			"usage : mpirun -n TOTAL_PROC %s mem -t 1 -o RESULTS REFERENCE FASTQ_R1 FASTQ_R2\n\n"
			"Requirements : After the creation of reference file with BWA you need to create a referenced\n"
			"	   map file genome with pidx like this  pidx ref.fasta generate a ref.fasta.map.\n"
			"	   The .map file is a copy of memory mapped reference used for shared memory purpose.\n",
			VERSION, progname);
		return 1; }

	/* Validate provided command (first argument) */
	if (strcmp(argv[1], "mem") != 0) {
		fprintf(stderr, "%s: unsupported %s command\n", progname, argv[1]);
		return 1; }

	/* initialize the BWA-MEM parameters to the default values */
	opt = mem_opt_init();
	memset(&opt0, 0, sizeof(opt0));
	while ((c = getopt(argc-1, argv+1, "1paMCSPVYjk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:N:W:x:G:h:y:K:X:H:o:")) >= 0) {
		if (c == 'k') opt->min_seed_len = atoi(optarg), opt0.min_seed_len = 1;
		else if (c == '1') ; /* FIXME: unsupported */
		else if (c == 'x') mode = optarg;
		else if (c == 'w') opt->w = atoi(optarg), opt0.w = 1;
		else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1;
		else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1;
		else if (c == 'O') {
			opt0.o_del = opt0.o_ins = 1;
			opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->o_ins = strtol(p+1, &p, 10);
		}
		else if (c == 'E') {
			opt0.e_del = opt0.e_ins = 1;
			opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->e_ins = strtol(p+1, &p, 10);
		}
		else if (c == 'L') {
			opt0.pen_clip5 = opt0.pen_clip3 = 1;
			opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->pen_clip3 = strtol(p+1, &p, 10);
		}
		else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1;
		else if (c == 'U') opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1;
		else if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		else if (c == 'p') opt->flag |= MEM_F_PE | MEM_F_SMARTPE;
		else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
		else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
		else if (c == 'Y') opt->flag |= MEM_F_SOFTCLIP;
		else if (c == 'V') opt->flag |= MEM_F_REF_HDR;
		else if (c == 'c') opt->max_occ = atoi(optarg), opt0.max_occ = 1;
		else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;
		else if (c == 'v') bwa_verbose = atoi(optarg);
		else if (c == 'j') ignore_alt = 1;
		else if (c == 'r') opt->split_factor = atof(optarg), opt0.split_factor = 1.;
		else if (c == 'D') opt->drop_ratio = atof(optarg), opt0.drop_ratio = 1.;
		else if (c == 'm') opt->max_matesw = atoi(optarg), opt0.max_matesw = 1;
		else if (c == 's') opt->split_width = atoi(optarg), opt0.split_width = 1;
		else if (c == 'G') opt->max_chain_gap = atoi(optarg), opt0.max_chain_gap = 1;
		else if (c == 'N') opt->max_chain_extend = atoi(optarg), opt0.max_chain_extend = 1;
		else if (c == 'W') opt->min_chain_weight = atoi(optarg), opt0.min_chain_weight = 1;
		else if (c == 'y') opt->max_mem_intv = atol(optarg), opt0.max_mem_intv = 1;
		else if (c == 'C') copy_comment = 1;
		else if (c == 'K') ; /* FIXME: unsupported */
		else if (c == 'X') opt->mask_level = atof(optarg);
		else if (c == 'h') {
			opt0.max_XA_hits = opt0.max_XA_hits_alt = 1;
			opt->max_XA_hits = opt->max_XA_hits_alt = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->max_XA_hits_alt = strtol(p+1, &p, 10);
		}
		else if (c == 'Q') {
			opt0.mapQ_coef_len = 1;
			opt->mapQ_coef_len = atoi(optarg);
			opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;
		}
		else if (c == 'R') {
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1;
		}
		else if (c == 'H') {
			if (optarg[0] != '@') {
				FILE *fp;
				if ((fp = fopen(optarg, "r")) != 0) {
					char *buf;
					buf = calloc(1, 0x10000);
					while (fgets(buf, 0xffff, fp)) {
						size_t i = strlen(buf);
						assert(buf[i-1] == '\n'); // a long line
						buf[i-1] = 0;
						hdr_line = bwa_insert_header(buf, hdr_line);
					}
					free(buf);
					fclose(fp);
				}
			} else hdr_line = bwa_insert_header(optarg, hdr_line);
		}
		else if (c == 'I') { // specify the insert size distribution
			pes0 = pes;
			pes[1].failed = 0;
			pes[1].avg = strtod(optarg, &p);
			pes[1].std = pes[1].avg * .1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].std = strtod(p+1, &p);
			pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);
			pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);
			if (pes[1].low < 1) pes[1].low = 1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].high = (int)(strtod(p+1, &p) + .499);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].low  = (int)(strtod(p+1, &p) + .499);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] mean insert size: %.3f, stddev: %.3f, max: %d, min: %d\n",
					__func__, pes[1].avg, pes[1].std, pes[1].high, pes[1].low);
		}
		/* Tool specific options */
		else if (c == 'o') file_out = optarg;
		else return 1; }

	if (mode) {
		if (strcmp(mode, "intractg") == 0) {
                        if (!opt0.o_del) opt->o_del = 16;
			if (!opt0.o_ins) opt->o_ins = 16;
			if (!opt0.b) opt->b = 9;
			if (!opt0.pen_clip5) opt->pen_clip5 = 5;
			if (!opt0.pen_clip3) opt->pen_clip3 = 5;
		} else if (strcmp(mode, "pacbio") == 0 || strcmp(mode, "pbref") == 0 || strcmp(mode, "ont2d") == 0) {
                        if (!opt0.o_del) opt->o_del = 1;
			if (!opt0.e_del) opt->e_del = 1;
			if (!opt0.o_ins) opt->o_ins = 1;
			if (!opt0.e_ins) opt->e_ins = 1;
			if (!opt0.b) opt->b = 1;
			if (opt0.split_factor == 0.) opt->split_factor = 10.;
			if (strcmp(mode, "ont2d") == 0) {
				if (!opt0.min_chain_weight) opt->min_chain_weight = 20;
				if (!opt0.min_seed_len) opt->min_seed_len = 14;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			} else {
				if (!opt0.min_chain_weight) opt->min_chain_weight = 40;
				if (!opt0.min_seed_len) opt->min_seed_len = 17;
				if (!opt0.pen_clip5) opt->pen_clip5 = 0;
				if (!opt0.pen_clip3) opt->pen_clip3 = 0;
			}
		} else {
			fprintf(stderr, "[E::%s] unknown read type '%s'\n", __func__, mode);
			return 1; // FIXME memory leak
		}
	} else {
		//update_a(opt, &opt0);
		if (opt0.a) { // matching score is changed
			if (!opt0.b) opt->b *= opt->a;
			if (!opt0.T) opt->T *= opt->a;
			if (!opt0.o_del) opt->o_del *= opt->a;
			if (!opt0.e_del) opt->e_del *= opt->a;
			if (!opt0.o_ins) opt->o_ins *= opt->a;
			if (!opt0.e_ins) opt->e_ins *= opt->a;
			if (!opt0.zdrop) opt->zdrop *= opt->a;
			if (!opt0.pen_clip5) opt->pen_clip5 *= opt->a;
			if (!opt0.pen_clip3) opt->pen_clip3 *= opt->a;
			if (!opt0.pen_unpaired) opt->pen_unpaired *= opt->a;
		}
	}

	nargs = argc - 1 - optind;
	if (opt->n_threads < 1) opt->n_threads = 1;
	if (nargs != 2 && nargs != 3) {
                fprintf(stderr, "\n");
		fprintf(stderr, "Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT        number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "       -k INT        minimum seed length [%d]\n", opt->min_seed_len);
		fprintf(stderr, "       -w INT        band width for banded alignment [%d]\n", opt->w);
		fprintf(stderr, "       -d INT        off-diagonal X-dropoff [%d]\n", opt->zdrop);
		fprintf(stderr, "       -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
		fprintf(stderr, "       -y INT        seed occurrence for the 3rd round seeding [%ld]\n", (long)opt->max_mem_intv);
		fprintf(stderr, "       -c INT        skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
		fprintf(stderr, "       -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [%.2f]\n", opt->drop_ratio);
		fprintf(stderr, "       -W INT        discard a chain if seeded bases shorter than INT [0]\n");
		fprintf(stderr, "       -m INT        perform at most INT rounds of mate rescues for each read [%d]\n", opt->max_matesw);
		fprintf(stderr, "       -S            skip mate rescue\n");
		fprintf(stderr, "       -P            skip pairing; mate rescue performed unless -S also in use\n");
		fprintf(stderr, "\nScoring options:\n\n");
		fprintf(stderr, "       -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [%d]\n", opt->a);
		fprintf(stderr, "       -B INT        penalty for a mismatch [%d]\n", opt->b);
		fprintf(stderr, "       -O INT[,INT]  gap open penalties for deletions and insertions [%d,%d]\n", opt->o_del, opt->o_ins);
		fprintf(stderr, "       -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [%d,%d]\n", opt->e_del, opt->e_ins);
		fprintf(stderr, "       -L INT[,INT]  penalty for 5'- and 3'-end clipping [%d,%d]\n", opt->pen_clip5, opt->pen_clip3);
		fprintf(stderr, "       -U INT        penalty for an unpaired read pair [%d]\n\n", opt->pen_unpaired);
		fprintf(stderr, "       -x STR        read type. Setting -x changes multiple parameters unless overriden [null]\n");
		fprintf(stderr, "                     pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)\n");
		fprintf(stderr, "                     ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)\n");
		fprintf(stderr, "                     intractg: -B9 -O16 -L5  (intra-species contigs to ref)\n");
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -p            smart pairing (ignoring in2.fq)\n");
		fprintf(stderr, "       -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n");
		fprintf(stderr, "       -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]\n");
		fprintf(stderr, "       -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "       -T INT        minimum score to output [%d]\n", opt->T);
		fprintf(stderr, "       -h INT[,INT]  if there are <INT hits with score >80%% of the max score, output all in XA [%d,%d]\n", opt->max_XA_hits, opt->max_XA_hits_alt);
		fprintf(stderr, "       -a            output all alignments for SE or unpaired PE\n");
		fprintf(stderr, "       -C            append FASTA/FASTQ comment to SAM output\n");
		fprintf(stderr, "       -V            output the reference FASTA header in the XR tag\n");
		fprintf(stderr, "       -Y            use soft clipping for supplementary alignments\n");
		fprintf(stderr, "       -M            mark shorter split hits as secondary\n\n");
		fprintf(stderr, "       -I FLOAT[,FLOAT[,INT[,INT]]]\n");
		fprintf(stderr, "                     specify the mean, standard deviation (10%% of the mean if absent), max\n");
		fprintf(stderr, "                     (4 sigma from the mean if absent) and min of the insert size distribution.\n");
		fprintf(stderr, "                     FR orientation only. [inferred]\n");
		fprintf(stderr, "\nExtra options:\n\n");
		fprintf(stderr, "       -o STR     the output file name\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "Note: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		free(opt);
		return 1;
	}
	bwa_fill_scmat(opt->a, opt->b, opt->mat);

	files = 0;
	file_ref = argv[optind+1+0];
	if (nargs > 1) {
		file_r1 = argv[optind+1+1]; files += 1;
	}
	if (nargs > 2) {
		file_r2 = argv[optind+1+2]; files += 1;
		opt->flag |= MEM_F_PE;
	}

	/* Check that output file (-o) is not null ... */
	if (file_out == NULL) {
		fprintf(stderr, "missing mandatory output file (-o)\n");
		res = MPI_Finalize();
		assert(res == MPI_SUCCESS);
		exit(2);
	}

	/* start up MPI */
	res = MPI_Init(&argc, &argv);
	assert(res == MPI_SUCCESS);
	res = MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	assert(res == MPI_SUCCESS);
	res = MPI_Comm_rank(MPI_COMM_WORLD, &rank_num);
	assert(res == MPI_SUCCESS);

	/* Derived file names */
	sprintf(file_map, "%s.map", file_ref);
	sprintf(file_tmp, "%s.offset.tmp", file_out);


	/* Work around build warning in non timing case */
	aft = 0; aft++;
	bef = 0; bef++;

	MPI_Get_processor_name(processor_name, &name_len);
	/*
	 * Map reference genome indexes in shared memory (by host)
	 */
	bef = MPI_Wtime();
	res = MPI_File_open(MPI_COMM_WORLD, file_map, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_map);
	assert(res == MPI_SUCCESS);
	res = MPI_File_get_size(fh_map, &size_map);
	assert(res == MPI_SUCCESS);
	res = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &comm_shr);
	assert(res == MPI_SUCCESS);
	res = MPI_Comm_rank(comm_shr, &rank_shr);
	assert(res == MPI_SUCCESS);
	size_shr = (rank_shr == 0) ? size_map : 0;
	res = MPI_Win_allocate_shared(size_shr, 0, MPI_INFO_NULL, comm_shr, &addr, &win_shr);
	assert(res == MPI_SUCCESS);
	res = MPI_Win_shared_query(win_shr, MPI_PROC_NULL, &size_shr, &res, &addr_map);
	assert(res == MPI_SUCCESS);

	//MPI_Win_lock_all(MPI_MODE_NOCHECK, win_shr);
	m = size_map; a = addr_map; size_tot = 0;
	while (rank_shr == 0) {
		res = MPI_File_read(fh_map, a, INT_MAX/2, MPI_UINT8_T, &status);
		assert(res == MPI_SUCCESS);
		res = MPI_Get_count(&status, MPI_UINT8_T, &count);
		assert(res == MPI_SUCCESS);
		if (count == 0) break;
		m -= count; a += count; size_tot += count;
	}
	assert(size_tot == 0 || size_tot == size_map);
	res = MPI_Win_fence(0, win_shr);
	assert(res == MPI_SUCCESS);

	bwa_mem2idx(size_map, addr_map, &indix);
	if (ignore_alt)
		for (c = 0; c < indix.bns->n_seqs; ++c)
			indix.bns->anns[c].is_alt = 0;

	res = MPI_File_close(&fh_map);
	assert(res == MPI_SUCCESS);
	aft = MPI_Wtime();
	xfprintf(stderr, "%s: mapped indexes (%.02f)\n", __func__, aft - bef);

	/*
	 * Create SAM header
	 * TODO: Add line for BWA version
	 */

	if (rank_num == 0) {
		int s, len;
		char *buff;
		MPI_File fh_out1;
		res = MPI_File_open(MPI_COMM_SELF, file_out, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_out1);
		assert(res == MPI_SUCCESS);
		res = MPI_File_set_size(fh_out1, 0);
		assert(res == MPI_SUCCESS);
		/* Add reference sequence lines */
		for (s = 0; s < indix.bns->n_seqs; ++s) {
			len = asprintf(&buff, "@SQ\tSN:%s\tLN:%d\n", indix.bns->anns[s].name, indix.bns->anns[s].len);
			res = MPI_File_write(fh_out1, buff, len, MPI_CHAR, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_Get_count(&status, MPI_CHAR, &count);
			assert(res == MPI_SUCCESS);
			assert(count == len);
			free(buff);
		}
		/* Add header lines */
		if (hdr_line != NULL) {
			len = asprintf(&buff, "%s\n", hdr_line);
			res = MPI_File_write(fh_out1, buff, len, MPI_CHAR, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_Get_count(&status, MPI_CHAR, &count);
			assert(res == MPI_SUCCESS);
			assert(count == len);
			free(buff);
		}
		/* Add read group line */
		if (rg_line != NULL) {
			len = asprintf(&buff, "%s\n", rg_line);
			res = MPI_File_write(fh_out1, buff, len, MPI_CHAR, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_Get_count(&status, MPI_CHAR, &count);
			assert(res == MPI_SUCCESS);
			assert(count == len);
			free(buff);
		}
		res = MPI_File_close(&fh_out1);
		assert(res == MPI_SUCCESS);
	}
	bef = MPI_Wtime();
	res = MPI_Barrier(MPI_COMM_WORLD);
	assert(res == MPI_SUCCESS);
	aft = MPI_Wtime();
	xfprintf(stderr, "%s: synched processes (%.02f)\n", __func__, aft - bef);

	rlen_r1 = rlen_r2 = READSIZE * opt->n_threads;

	seqs = NULL;
	uint64_t bases = opt->chunk_size * opt->n_threads;
	curr[0] =  curr[1] = curr[2] = curr[3] =  curr[4] = curr[5] = 0;
	int init;
	off_t start_off_r1 = 0;
	off_t start_off_r2 = 0;

	/* Open temporary offset file in appending mode*/
	res = MPI_File_open(MPI_COMM_WORLD, file_tmp, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fd_tmp);
	assert(res == MPI_SUCCESS);

	/* initialize temporary offset file */
	if (rank_num == 0){

		res = MPI_File_write_at(fd_tmp, 0, curr, 6, MPI_LONG_LONG_INT, &status);
		assert(res == MPI_SUCCESS);
		res = MPI_File_sync(fd_tmp) ;
		assert(res == MPI_SUCCESS);
	}

	/* Open input files */
	res = MPI_File_open(MPI_COMM_WORLD, file_r1, MPI_MODE_RDWR, MPI_INFO_NULL, &fh_r1);
	assert(res == MPI_SUCCESS);

	res = MPI_File_open(MPI_COMM_WORLD, file_r2, MPI_MODE_RDWR, MPI_INFO_NULL, &fh_r2);
	assert(res == MPI_SUCCESS);

	res = MPI_File_open(MPI_COMM_WORLD, file_r1, MPI_MODE_RDWR, MPI_INFO_NULL, &fh_r3);
	assert(res == MPI_SUCCESS);

	res = MPI_File_open(MPI_COMM_WORLD, file_r2, MPI_MODE_RDWR, MPI_INFO_NULL, &fh_r4);
	assert(res == MPI_SUCCESS);

	/* Open output file */
	res = MPI_File_open(MPI_COMM_WORLD, file_out, MPI_MODE_WRONLY|MPI_MODE_APPEND, MPI_INFO_NULL, &fh_out);
	assert(res == MPI_SUCCESS);

	rma_mutex_init_offset(&mutex_offset, MPI_COMM_WORLD, 0);

	MPI_Barrier(MPI_COMM_WORLD);

	size_t total_line_read_offset_proc = 0;
	size_t total_line_read_align_proc = 0;

	if (rank_num < NUM_PROC_OFFSET){

		int *tmp_tab = malloc(proc_num * sizeof(int));
		int k1 = 0;

		for (k1 = 0; k1 < proc_num; k1++)
		tmp_tab[k1] = k1%NUM_PROC_OFFSET;
		int num_consumer = 0;

		for (k1 = NUM_PROC_OFFSET; k1 < proc_num; k1++){
			if (tmp_tab[k1] == rank_num) num_consumer++;
		}

		fprintf(stderr, "Rank %d :: number of consumer = % d \n", rank_num, num_consumer);
		int consumer = 0;

		while (1) {

			rma_mutex_lock_offset(&mutex_offset);

			char *x, *b1, *b2;
			size_t csiz = sizeof(curr);
			ssize_t ssiz;
			size_t size_r1, size_r2;
			int count_r1, full_r1, line_number;
			uint64_t seqs_r1, base_r1, seqs_r2, base_r2;
			ptrdiff_t diff_r1, diff_r2;
			bef = MPI_Wtime();

			res = MPI_File_read_at(fd_tmp, 0, curr, 6, MPI_LONG_LONG_INT, &status);
			assert(res == MPI_SUCCESS);

			start_off_r1 = curr[0] + curr[1];
			start_off_r2 = curr[3] + curr[4];

			/* Read current chunk data to be processed ... */
	again:

			if (file_r1 != NULL) {

				if (buffer_r1) buffer_r1 = (char *)realloc(buffer_r1, rlen_r1 + 1);
				else buffer_r1 = malloc(rlen_r1+1);
				buffer_r1[rlen_r1] = '\0';
				assert(buffer_r1 != NULL);
				res = MPI_File_read_at(fh_r1, start_off_r1, buffer_r1, (int)rlen_r1, MPI_CHAR, &status);
				assert(res == MPI_SUCCESS);
				res = MPI_Get_count(&status, MPI_CHAR, &count_r1);
				assert(res == MPI_SUCCESS);
				assert(*buffer_r1 == '@');
			}

			if (file_r2 != NULL) {

				if (buffer_r2) buffer_r2 = (char *)realloc(buffer_r2, rlen_r2 + 1);
				else buffer_r2 = malloc(rlen_r2+1);
				buffer_r2[rlen_r2] = '\0';
				assert(buffer_r2 != NULL);
				res = MPI_File_read_at(fh_r2, start_off_r2, buffer_r2, (int)rlen_r2, MPI_CHAR, &status);
				assert(res == MPI_SUCCESS);
				res = MPI_Get_count(&status, MPI_CHAR, &count_r2);
				assert(res == MPI_SUCCESS);
				assert(*buffer_r2 == '@');
			}

			/* Check for expected base count */
			b1 = buffer_r1; size_r1 = (size_t)count_r1; line_r1 = 0; base_r1 = 0;
			b2 = buffer_r2; size_r2 = (size_t)count_r2; line_r2 = 0; base_r2 = 0;
			while (1) {
				off_t val1 = 0;
				off_t val2 = 0;
				int zero = 0;
				if (file_r1 != NULL) {
					x = memchr(b1, '\n', size_r1);
					if (x == NULL) break;
					diff_r1 = x - b1 + 1;
					b1 = x + 1;
					line_r1++;
					size_r1 -= diff_r1;

					if (line_r1 % 4 == 0) {
						base_r1 += (uint64_t)diff_r1 - 1;
					}
				}
				if (file_r2 != NULL) {
					x = memchr(b2, '\n', size_r2);
					if (x == NULL) break;
					diff_r2 = x - b2 + 1; b2 = x + 1;
					line_r2++; size_r2 -= diff_r2;
					if (line_r2 % 4 == 0) {
						base_r2 += (uint64_t)diff_r2 - 1; }
					}
				if (base_r1 + base_r2 >= bases) break;
			}

			assert(line_r2 == line_r1);
			full_r1 = (file_r1 != NULL) ? count_r1 == rlen_r1 : 0;
			full_r2 = (file_r2 != NULL) ? count_r2 == rlen_r2 : 0;

			if ((full_r1 || full_r2) && ((base_r1 + base_r2) < bases)) {
					rlen_r1 += READSIZE * opt->n_threads;
					rlen_r2 += READSIZE * opt->n_threads;
					goto again;
			}
			size_r1 = b1 - buffer_r1;
			assert(size_r1 <= rlen_r1);
			size_r2 = b2 - buffer_r2;
			assert(size_r2 <= rlen_r2);
			aft = MPI_Wtime();

			assert(size_r1 == size_r2);


			total_line_read_offset_proc += line_r2;
			total_line_read_offset_proc += line_r1;

			/* Stop if nothing to process */
			if (size_r1 == 0 || size_r2 == 0){
				rma_mutex_unlock_offset(&mutex_offset);
				break;
			}

			curr[0] = start_off_r1;
			curr[1] = size_r1;
			curr[2] = line_r1;
			curr[3] = start_off_r2;
			curr[4] = size_r2;
			curr[5] = line_r2;

			res = MPI_File_write_at(fd_tmp, 0, curr, 6, MPI_LONG_LONG_INT, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_File_sync(fd_tmp) ;
			assert(res == MPI_SUCCESS);

			rma_mutex_unlock_offset(&mutex_offset);

			//we recieve the consumer rank
			res = MPI_Recv(&consumer, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			assert(res == MPI_SUCCESS);
			MPI_Datatype send_buf;
			res = MPI_Type_contiguous (6, MPI_LONG_LONG_INT, &send_buf);
			assert(res == MPI_SUCCESS);
			res = MPI_Type_commit(&send_buf);
			assert(res == MPI_SUCCESS);
			res = MPI_Send(&curr , 1, send_buf, consumer, 0, MPI_COMM_WORLD);
			assert(res == MPI_SUCCESS);
			res = MPI_Type_free( &send_buf);
			assert(res == MPI_SUCCESS);
			//should not come here!!
			if ((file_r1 != NULL && size_r1 == 0))
				fprintf(stderr, "[M::%s] sequences count differ between R1 and R2\n", __func__);

		}//end while(1)

		//here we shall tell the jobs we are finished

		int t=0;
		size_t curr_zero[6];
		curr_zero[0] = 0;
		curr_zero[1] = 0;
		curr_zero[2] = 0;
		curr_zero[3] = 0;
		curr_zero[4] = 0;
		curr_zero[5] = 0;

		do{

			res = MPI_Recv(&consumer, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			assert(res == MPI_SUCCESS);
			//then we send to the consumer our offsets
			MPI_Datatype send_buf;
			res = MPI_Type_contiguous (6, MPI_LONG_LONG_INT, &send_buf);
			assert(res == MPI_SUCCESS);
			res = MPI_Type_commit(&send_buf);
			assert(res == MPI_SUCCESS);
			res = MPI_Send(&curr_zero , 1, send_buf, consumer, 0, MPI_COMM_WORLD);
			assert(res == MPI_SUCCESS);
			res = MPI_Type_free(&send_buf);
			assert(res == MPI_SUCCESS);

			t++;
		}while (t < num_consumer);
	}

	else {

		int prov = 0;
		prov = rank_num % NUM_PROC_OFFSET;
		int write_flag = 0;
		size_t buffer_out_size = 0;

		while (1) {

			res = MPI_Send(&rank_num , 1, MPI_INT, prov, 0, MPI_COMM_WORLD);
			assert(res == MPI_SUCCESS);

			MPI_Datatype recv_buf;
			res = MPI_Type_contiguous (6, MPI_LONG_LONG_INT, &recv_buf);
			assert(res == MPI_SUCCESS);
			res = MPI_Type_commit(&recv_buf);
			assert(res == MPI_SUCCESS);
			//we recieve the provider offsets
			res = MPI_Recv(&curr, 1, recv_buf, prov, 0, MPI_COMM_WORLD, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_Get_count( &status, recv_buf, &count );
			assert(res == MPI_SUCCESS);
			assert(count == 1);
			res = MPI_Type_free( &recv_buf );
			assert(res == MPI_SUCCESS);

			if (curr[4] == 0 || curr[1] == 0 || curr[2] == 0 || curr[5] == 0) break; //no line to read it's finished
			bef = MPI_Wtime();

			if (file_r1 != NULL) {
				/* Read some sequences */

				if (buffer_r1) buffer_r1 = (char *)realloc(buffer_r1, curr[1]+1);
				else buffer_r1 = malloc(curr[1]+1);
				assert(buffer_r1 != NULL);
				buffer_r1[curr[1]] = '\0';

				res = MPI_File_read_at(fh_r3, (MPI_Offset)curr[0], buffer_r1, (int)curr[1], MPI_CHAR, &status);
				assert(res == MPI_SUCCESS);
				res = MPI_Get_count(&status, MPI_CHAR, &count);
				size_r1 = (size_t)count;
				assert(res == MPI_SUCCESS);
				assert(*buffer_r1 == '@');
			}
			if (file_r2 != NULL) {

				if (buffer_r2) buffer_r2 = (char *)realloc(buffer_r2, curr[4]+1);
				else buffer_r2 = malloc(curr[4]+1);
				assert(buffer_r2 != NULL);
				buffer_r2[curr[4]] = '\0';

				res = MPI_File_read_at(fh_r4, (MPI_Offset)curr[3], buffer_r2, (int)curr[4], MPI_CHAR, &status);
				assert(res == MPI_SUCCESS);
				res = MPI_Get_count(&status, MPI_CHAR, &count);
				size_r2 = (size_t)count;
				assert(res == MPI_SUCCESS);
				assert(*buffer_r2 == '@');
			}

			aft = MPI_Wtime();
			fprintf(stderr, "Rank %d: read size of %zu file (%.02f)\n", rank_num, strlen(buffer_r2),aft - bef);

			size_t line_r1 = curr[2];
			size_t line_r2 = curr[5];
			/* Allocate sequence structure */
			assert(line_r1 % 4 == 0);
			assert(line_r2 % 4 == 0);
			size_t seqs_r1 = line_r1 / 4;
			size_t seqs_r2 = line_r2 / 4;
			seqs = malloc((seqs_r1 + seqs_r2) * sizeof(*seqs));
			assert(seqs != NULL);

			/* Parse sequences ... */
			bef = MPI_Wtime();
			if (file_r1 != NULL) {
				p = q = buffer_r1; e = buffer_r1 + size_r1; line_number = 0;
				while (q < e) {
					if (*q != '\n') { q++; continue; }
					/* We have a full line ... process it */
					*q = '\0';
					n = files * (line_number / 4);
					switch (line_number % 4) {
					case 0: /* Line1: Name and Comment */
						assert(*p == '@');
						seqs[n].name = p + 1;
						while (*p && !isspace((unsigned char)*p)) p++;
						if (*(p-2) == '/' && isdigit((unsigned char)*(p-1))) *(p-2) = '\0';
						if (*p) *p++ = '\0';
						seqs[n].comment = (copy_comment != 0) ? p : NULL;
						seqs[n].sam = NULL;
						break;
					case 1: /* Line2: Sequence */
						seqs[n].seq = p;
						seqs[n].l_seq = q - p;
						break;
					case 2: /* Line3: Ignored */
						assert(*p == '+');
						break;
					case 3: /* Line4: Quality */
						seqs[n].qual = p;
						break; }
					p = ++q; line_number++;
				}
				assert(curr[2] == line_number);
			}


			if (file_r2 != NULL) {
				p = q = buffer_r2; e = buffer_r2 + size_r2; line_number = 0;
				while (q < e) {
					if (*q != '\n') { q++; continue; }
					/* We have a full line ... process it */
					*q = '\0'; n = files * (line_number / 4) + 1;
					switch (line_number % 4) {
					case 0: /* Line1: Name and Comment */
						assert(*p == '@');
						seqs[n].name = p + 1;
						while (*p && !isspace((unsigned char)*p)) p++;
						if (*(p-2) == '/' && isdigit((unsigned char)*(p-1))) *(p-2) = '\0';
						if (*p) *p++ = '\0';
						seqs[n].comment = (copy_comment != 0) ? p : NULL;
						seqs[n].sam = NULL;
						break;
					case 1: /* Line2: Sequence */
						seqs[n].seq = p;
						seqs[n].l_seq = q - p;
						break;
					case 2: /* Line3: Ignored */
						assert(*p == '+');
						break;
					case 3: /* Line4: Quality */
						seqs[n].qual = p;
						break; }
					p = ++q; line_number++;
				}
				assert(curr[5] == line_number);
			}

			total_line_read_align_proc += curr[5];
			total_line_read_align_proc += curr[2];

			aft = MPI_Wtime();

			/* Process sequences */
			bef = MPI_Wtime();
			mem_process_seqs(opt, indix.bwt, indix.bns, indix.pac, 0, (int)(seqs_r1+seqs_r2), seqs, pes0);
			aft = MPI_Wtime();
			fprintf(stderr, "Rank %d: computed mappings (%.02f)\n", rank_num, aft - bef);
			/* FIXME: Write results */
			bef = MPI_Wtime();
			localsize = 0;

			for (n = 0; n < seqs_r1+seqs_r2; n++) {
				/* Reuse .l_seq to store SAM line length to avoid multiple strlen() calls */
				seqs[n].l_seq = strlen(seqs[n].sam);
				localsize += seqs[n].l_seq;
			}

			assert(localsize <= INT_MAX);
			buffer_out = (char *)realloc(buffer_out,  localsize +1);
			assert(buffer_out != NULL);
			buffer_out[localsize] = '\0';
			p = buffer_out;

			for (n = 0; n < seqs_r1+seqs_r2; n++) {
				memmove(p, seqs[n].sam, seqs[n].l_seq);
				p += seqs[n].l_seq;
				free(seqs[n].sam);
			}

			res = MPI_File_write_shared(fh_out, buffer_out, (int)localsize, MPI_CHAR, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_Get_count(&status, MPI_CHAR, &count);
			assert(res == MPI_SUCCESS);
			assert(count == (int)localsize);

			free(seqs);

		}
		fprintf(stderr, "Rank %d: finish aligning \n", rank_num);
	}

	res = MPI_Barrier(MPI_COMM_WORLD);
	assert(res == MPI_SUCCESS);


	/*
	 * For DEBUG
	 *
	 * Get the number of line read and bytes writen
	 *
	MPI_Allreduce(&total_writen_by_proc, &total_writen, proc_num, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
	fprintf(stderr, "Rank %d: total writen by proc= %zu \n", rank_num, total_writen_by_proc);
	if (total_line_read_offset_proc)
		fprintf(stderr, "Rank %d: total line read proc offset = %zu \n", rank_num, total_line_read_offset_proc);

	if (total_line_read_align_proc)
		fprintf(stderr, "Rank %d: total line read proc align = %zu \n", rank_num,  total_line_read_align_proc);
	*/

	rma_mutex_free_offset(&mutex_offset);

	res = MPI_File_close(&fh_r1);
	assert(res == MPI_SUCCESS);
	res = MPI_File_close(&fh_r2);
	assert(res == MPI_SUCCESS);
	res = MPI_File_close(&fd_tmp);
	assert(res == MPI_SUCCESS);
	res = MPI_File_close(&fh_r4);
	assert(res == MPI_SUCCESS);
	res = MPI_File_close(&fh_r3);
	assert(res == MPI_SUCCESS);
	res = MPI_File_close(&fh_out);
	assert(res == MPI_SUCCESS);

	if (buffer_out) free(buffer_out);
	if (buffer_r2)  free(buffer_r2);
	if (buffer_r1)  free(buffer_r1);

	if (opt) free(opt);

	(void)unlink(file_tmp);

	res = MPI_Win_free(&win_shr);
	assert(res == MPI_SUCCESS);
	res = MPI_Comm_free(&comm_shr);
	assert(res == MPI_SUCCESS);

	res = MPI_Barrier(MPI_COMM_WORLD);
	assert(res == MPI_SUCCESS);
	fprintf(stderr, "rank %d we are finished \n", rank_num);
	res = MPI_Finalize();
	assert(res == MPI_SUCCESS);
	return 0;
}

