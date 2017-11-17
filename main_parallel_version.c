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

#include <sys/mman.h>
#include <sys/stat.h>

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

#ifdef TIMIMG
#define xfprintf fprintf
#else
#define xfprintf(...) /**/
#endif



void init_goff(size_t *goff, MPI_File mpi_filed, size_t fsize,int numproc,int rank){


	char * current_line = NULL;
	MPI_Status status;
	int i = 0;
	int j = 0;
	// TODO problem here when we have only one job!
	size_t lsize = fsize/numproc;
	goff[0]=0;
	for(i=1;i<numproc;i++){goff[i]=lsize*i;}
	goff[numproc]=fsize;
	for(i=1;i<numproc;i++){
		current_line =(char*)calloc(2000,sizeof(char));
		//we read 2K caracters
		MPI_File_read_at(mpi_filed, (MPI_Offset)goff[i], current_line, 2000, MPI_CHAR, &status);
		j=0;
		while (j < 2000 ){
			//we check the \n after and before
			j++;
			if (current_line[j] == '+' && current_line[j+1] == '\n' && current_line[j-1] == '\n') {j++;break;}}
		//then we look for return after quality
		while( j < 2000){j++;if (current_line[j] == '\n') {break;}}
		goff[i]+=(j+1); free(current_line);}

	return;
}

int main(int argc, char *argv[]) {
	double bef_total_run,aft_total_run,bef, aft;

	char *file_r1 = NULL, *file_r2 = NULL;
	struct stat stat_r1, stat_r2;
	char *buffer_r1, *buffer_r2;
	int fd_in1, fd_in2;
	uint8_t *a, *addr, *addr_map;
	char *file_out = NULL;
	char *buffer_out;
	char *file_ref = NULL;
	char file_map[PATH_MAX], file_tmp[PATH_MAX];
	int proc_num, rank_num, rank_shr;
	size_t localsize;
	size_t n = 0;

	MPI_Aint size_shr;
	MPI_Comm comm_shr;
	MPI_File fh_map, fh_r1, fh_r2, fh_out, fh_tmp, fh_tmp1;
	MPI_Offset *coff, m, size_map, size_tot;
	MPI_Status status;
	MPI_Win win_shr;

	mem_opt_t *opt, opt0;
	mem_pestat_t pes[4], *pes0 = NULL;
	bwaidx_t indix;
	bseq1_t *seqs;

	int res, count;
	int files, nargs;

	char *p, *q, *s, *e;
	off_t locoff, globoff, locsiz, *alloff, *curoff, maxsiz, chunck_sz, totsiz, filsiz;

	int c, copy_comment = 0;
	char *rg_line = NULL, *hdr_line = NULL;
	int ignore_alt = 0;
	const char *mode = NULL;

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

	/* Derived file names */
	sprintf(file_map, "%s.map", file_ref);
	sprintf(file_tmp, "%s.tmp", file_out);

	/* start up MPI */
	res = MPI_Init(&argc, &argv);
	assert(res == MPI_SUCCESS);
	res = MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	assert(res == MPI_SUCCESS);
	res = MPI_Comm_rank(MPI_COMM_WORLD, &rank_num);
	assert(res == MPI_SUCCESS);

	/* Check that output file (-o) is not null ... */
	if (file_out == NULL) {
		fprintf(stderr, "missing mandatory output file (-o)\n");
		res = MPI_Finalize();
		assert(res == MPI_SUCCESS);
		exit(2);
	}

	/* Check R1 & R2 file sizes */
	if (file_r1 != NULL && stat(file_r1, &stat_r1) == -1) {
		fprintf(stderr, "%s: %s\n", file_r1, strerror(errno));
		res = MPI_Finalize();
		assert(res == MPI_SUCCESS);
		exit(2);
	}
	if (file_r2 != NULL && stat(file_r2, &stat_r2) == -1) {
		fprintf(stderr, "%s: %s\n", file_r2, strerror(errno));
		res = MPI_Finalize();
		assert(res == MPI_SUCCESS);
		exit(2);
	}
	if (file_r2 != NULL && stat_r1.st_size != stat_r2.st_size) {
		fprintf(stderr, "file sizes for R1 and R2 do not match\n");
		//res = MPI_Finalize();
		//assert(res == MPI_SUCCESS);
		//exit(2);
	}

	/* Work around build warning in non timing case */
	aft = 0; aft++;
	bef = 0; bef++;
	bef_total_run = 0; bef_total_run++;
	aft_total_run  = 0; aft_total_run++;
	/*
	 * Map reference genome indexes in shared memory (by host)
	 */
	bef_total_run = MPI_Wtime();
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

	m = size_map; a = addr_map; size_tot = 0;
	while (rank_shr == 0) {
	    res = MPI_File_read(fh_map, a, INT_MAX/2, MPI_UINT8_T, &status);
	    assert(res == MPI_SUCCESS);
	    res = MPI_Get_count(&status, MPI_UINT8_T, &count);
	    assert(res == MPI_SUCCESS);
	    if (count == 0) break;
	    m -= count; a += count; size_tot += count; }
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

	    res = MPI_File_delete(file_out, MPI_INFO_NULL);
	    assert(res == MPI_SUCCESS || res == MPI_ERR_NO_SUCH_FILE || res == MPI_ERR_IO);
	    res = MPI_File_open(MPI_COMM_SELF, file_out, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_out);
	    assert(res == MPI_SUCCESS);
	    /* Add reference sequence lines */
	    for (s = 0; s < indix.bns->n_seqs; ++s) {
		len = asprintf(&buff, "@SQ\tSN:%s\tLN:%d\n", indix.bns->anns[s].name, indix.bns->anns[s].len);
		res = MPI_File_write(fh_out, buff, len, MPI_CHAR, &status);
		assert(res == MPI_SUCCESS);
		res = MPI_Get_count(&status, MPI_CHAR, &count);
		assert(res == MPI_SUCCESS);
		assert(count == len);
		free(buff);
	    }
	    /* Add header lines */
	    if (hdr_line != NULL) {
		len = asprintf(&buff, "%s\n", hdr_line);
		res = MPI_File_write(fh_out, buff, len, MPI_CHAR, &status);
		assert(res == MPI_SUCCESS);
		res = MPI_Get_count(&status, MPI_CHAR, &count);
		assert(res == MPI_SUCCESS);
		assert(count == len);
		free(buff);
	    }
	    /* Add read group line */
	    if (rg_line != NULL) {
		len = asprintf(&buff, "%s\n", rg_line);
		res = MPI_File_write(fh_out, buff, len, MPI_CHAR, &status);
		assert(res == MPI_SUCCESS);
		res = MPI_Get_count(&status, MPI_CHAR, &count);
		assert(res == MPI_SUCCESS);
		assert(count == len);
		free(buff);
	    }
	    res = MPI_File_close(&fh_out);
	    assert(res == MPI_SUCCESS);
	}
	bef = MPI_Wtime();
	res = MPI_Barrier(MPI_COMM_WORLD);
	assert(res == MPI_SUCCESS);
	aft = MPI_Wtime();
	xfprintf(stderr, "%s: synched processes (%.02f)\n", __func__, aft - bef);

	/******* Now we split in 2 cases ************/

	if (stat_r1.st_size == stat_r2.st_size) {
		
		/*
		 * We are in the case the reads are not trimmed
		 *  First we divide the file in big windows. Each window is a size of 1Gb * number of jobs 
		 * Then on each window we split again per jobs. 
		 * Finally each jobs compute offsets of chunks of 1Mb approximately 
		 * /

		
		/* The rank 0 estimates the chunck size*/
		off_t tmp_sz     = 1024; 
		if (rank_num == 0){		
			// 512 Mo are  
			int fd_tmp = open(file_r1, O_RDONLY, 0666);
			char *buffer = malloc(tmp_sz  + 1);
			buffer[tmp_sz] = '\0';	
			read(fd_tmp, buffer, tmp_sz);
			assert(strlen(buffer) == tmp_sz);
			assert( *buffer == '@');
			/* Estimate sequence size (bases and bytes) */
			/* With this estimation we compute approximately the chuncks size and offset*/ 
			size_t slen, blen;
			s = buffer;
			e = buffer + tmp_sz; 
			p = q = s;
			while (q < e && *q != '\n') q++; p = ++q;
			while (q < e && *q != '\n') q++; blen = q - p; p = ++q;
			while (q < e && *q != '\n') q++; p = ++q;
			while (q < e && *q != '\n') q++; slen = q - s + 1;
		
			/* Split local buffer in chunks of 100000000 bases */
			chunck_sz = 80000000 / blen * slen;
			//chunck_sz = 1000000 / blen * slen;
			chunck_sz *= opt->n_threads;
			chunck_sz /= files; 		
			assert(close(fd_tmp) != -1);
			free(buffer);		
		}
		
		/* the rank 0 broadcast the */
		res = MPI_Bcast(&chunck_sz, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD); 
		assert(res == MPI_SUCCESS);
		//fprintf(stderr, "%s: rank %d : chunck_sz = %zu \n", __func__, rank_num, chunck_sz);
		/*we divide the file in window of chunk_sz * nb_job*/
		off_t file_sz    = stat_r1.st_size;
		size_t total_local_reads = 0;     
		size_t total_reads = 0;
		off_t total_sz    = 0;
		int k = 0;
		off_t *window_offsets = NULL;
		off_t end_offset_window = 0;
		int number_of_window = 0;
		
		// Each job job compute offsets of the chuncks of chunck_sz buffer
		off_t window_sz  =  chunck_sz * proc_num; 		//total window size is chunck_size * number of jobs
		if ( window_sz > file_sz){
			number_of_window = 0;
			window_offsets = calloc( 2, sizeof(off_t));
			window_sz = file_sz;
			end_offset_window = file_sz;
			window_offsets[1] = file_sz;	
		}
		else  {
			number_of_window = stat_r1.st_size / window_sz;
			fprintf(stderr, "%s: rank %d : number of windows = %zu \n", __func__, rank_num, number_of_window);
			window_offsets = calloc(number_of_window+ 2, sizeof(off_t)); 
			end_offset_window = (rank_num + 1) * window_sz;
		}				
		  
		//Rank 0 compute the begin and stop of each window
		// TODO: optimize this part involving more jobs
		if (number_of_window > 0){

			int fd_in0 = open(file_r1, O_RDONLY, 0666);
			char *buffer_r0 = malloc( tmp_sz + 1);
			buffer_r0[tmp_sz] = '\0';	
			if ( rank_num == 0 ){
			
				for (k = 0; k < (number_of_window + 1); k++){ window_offsets[k] = window_sz * k;}
				window_offsets[number_of_window + 1] = file_sz;			
				/* now we adjust the offset of the window to start at the beginning of a read*/ 
				for (k = 1; k < (number_of_window + 1); k++){
					end_offset_window = window_offsets[k];
					pread(fd_in0, buffer_r0, tmp_sz, end_offset_window);
					assert(strlen(buffer_r0) == tmp_sz);
					p = buffer_r0;
					e = buffer_r0 + tmp_sz;
					while (p < e) {
						if (*p != '@') { p++; continue; }
						if (p != buffer_r0 && *(p-1) != '\n') { p++; continue; }
						q = p + 1;
						while (q < e && *q != '\n') q++; q++;
						while (q < e && *q != '\n') q++; q++;
						if (q < e && *q == '+') break;
						p++;
					}
					assert( *p == '@');
					//we update begining offsets of windows in order to start at the begining of a read			
					 window_offsets[k]   += p - buffer_r0; 
				}
			}
			// rank0 scatter the result vector of windows offsets		
			res = MPI_Bcast( window_offsets, (number_of_window + 2), MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD); 
			assert(res == MPI_SUCCESS);
			assert(close(fd_in0) != -1);
			free(buffer_r0);		
		}
		
		 /* In each window we compute the begining offsets of the chuncks and we save it */
		off_t tmp_off = 0; 
		total_sz  = 0;
		fd_in1 = open(file_r1, O_RDONLY, 0666);
		assert(fd_in1 != -1);
		res = MPI_File_delete(file_tmp, MPI_INFO_NULL);
		assert(res == MPI_SUCCESS || res == MPI_ERR_NO_SUCH_FILE || res == MPI_ERR_IO);
		res = MPI_File_open(MPI_COMM_WORLD, file_tmp, MPI_MODE_WRONLY| MPI_MODE_CREATE| MPI_MODE_APPEND, MPI_INFO_NULL, &fh_tmp);
		assert(res == MPI_SUCCESS);
		size_t total_chuncks = 0;
		size_t total_window_sz = 0;
		alloff = calloc((size_t)proc_num + 1, sizeof(*alloff));
		assert(alloff != NULL);

		for ( k = 0; k < (number_of_window + 1); k ++){	
			
			/*
			First we search the begining of the read in a local buffer
			each rank read 500Mo of fastq and go to the first read QName
			*/
			bef = MPI_Wtime();
			off_t window_sz  =window_offsets[k+1] - window_offsets[k];
			
			//locsiz = window_sz / proc_num ;
			locsiz = chunck_sz;
			globoff = rank_num * locsiz + window_offsets[k];		
			
			if (globoff > window_offsets[k+1]) {
				/* the last chunck could be outside of the file*/
				globoff =  window_offsets[k+1];
				tmp_sz = 0;		
			}
			else{
				/* otherwise the */ 
				buffer_r1 = malloc(tmp_sz + 1);
				buffer_r1[tmp_sz] = '\0';	
				/*buffer_r1 contains 500Mo of reads fastq*/
				pread(fd_in1, buffer_r1, tmp_sz, globoff);
				assert(strlen(buffer_r1) == tmp_sz);	
				/*we look for the first QNAME read*/
				if ( rank_num == 0 ) assert( *buffer_r1 == '@'); //make sure the window start on a QName
				p = buffer_r1;
				e = buffer_r1 + tmp_sz;
				while (p < e) {
					if (*p != '@') { p++; continue; }
					if (p != buffer_r1 && *(p-1) != '\n') { p++; continue; }
					q = p + 1;
					while (q < e && *q != '\n') q++; q++;
					while (q < e && *q != '\n') q++; q++;
					if (q < e && *q == '+') break;
					p++;
					}
				assert( *p == '@');
				//we update the global offset in the fastq file			
				globoff  += p - buffer_r1;
			}			
			//we save the new global offset in a vector alloff
			curoff = alloff + proc_num; 
			*curoff =window_offsets[k+1] ;
			res = MPI_Allgather(&globoff, 1, MPI_LONG_LONG_INT, alloff, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
			assert(res == MPI_SUCCESS);
			curoff = alloff + rank_num; 
			locsiz = *(curoff+1) - *(curoff);
			
			/* in case a rank has nothing to process it writes nothing*/
			if ( locsiz != 0){ 
				/*we free buffer_r1 and reload it*/						
				size_t total_chucks_sz = locsiz;
				buffer_r1 = realloc(buffer_r1, locsiz + 1);
				buffer_r1[locsiz] = '\0';	
				pread(fd_in1, buffer_r1, locsiz, *curoff);
				assert(strlen(buffer_r1) == locsiz);
				assert( *buffer_r1 == '@');
						
				totsiz = 0;
				filsiz = locsiz / chunck_sz + 1;
				res = MPI_Type_size(MPI_OFFSET, &c);
				assert(res == MPI_SUCCESS);
				coff = malloc(filsiz * 2 * c);
				assert(coff != NULL);
				c = 0;
				locoff = *curoff;
				s = buffer_r1;
				e = buffer_r1 + locsiz;
				maxsiz = chunck_sz; 
				while (1) {
					if (e - s < maxsiz) { maxsiz = e - s; }
					p = s + maxsiz;
					while (p < e) {
						if (*p != '@') { p++; continue; }
						if (p != buffer_r1 && *(p-1) != '\n') { p++; continue; }
						q = p + 1;
						while (q < e && *q != '\n') q++; q++;
						while (q < e && *q != '\n') q++; q++;
						if (q < e && *q == '+') break;
						p++; 
					}
					locsiz = p - s;
					if (locsiz == 0) break;
					coff[c++] = locoff; 
					coff[c++] = locsiz;
					s += locsiz; 
					locoff += locsiz; 
					totsiz += locsiz;
				 }
				assert (total_chucks_sz == totsiz);
			}
			else{
				
				filsiz = 0;
				totsiz = 0;
				res = MPI_Type_size(MPI_OFFSET, &c);
				assert(res == MPI_SUCCESS);
				coff = malloc(2 * c);
				assert(coff !=NULL );
				c = 0;
				coff[c++] = 0;
				coff[c++] = 0;	
				//c++;								
			}
			MPI_Barrier(MPI_COMM_WORLD);
			total_chuncks += c;
			/* Save offsets+sizes in temp file */
			
			res = MPI_File_write_ordered(fh_tmp, coff, c, MPI_OFFSET, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_Get_count(&status, MPI_OFFSET, &count);
			assert(res == MPI_SUCCESS);
			assert(count == c);
			free(coff);
			if (locsiz != 0) free(buffer_r1);	
			
			res = MPI_Allreduce(&totsiz, &filsiz, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
			assert(res == MPI_SUCCESS);			
			total_window_sz += window_sz;	
			aft = MPI_Wtime();
			fprintf(stderr, "%s: rank %d : Time to compute chunks (%.02f)\n", __func__, rank_num, aft - bef);
			assert(filsiz == window_sz);
	
		} 

		free(alloff); 
		// end loop computing chuncks		
		assert(total_window_sz == stat_r1.st_size);
		assert(close(fd_in1) != -1);
		fprintf(stderr, "%s: rank %d : total chunks = %zu \n", __func__, rank_num, total_chuncks);
		
		MPI_Offset file_offset_sz = 0;
		res = MPI_File_get_size(fh_tmp, &file_offset_sz);
		assert(res == MPI_SUCCESS);
		MPI_Barrier(MPI_COMM_WORLD);

		res = MPI_File_close(&fh_tmp);
		assert(res == MPI_SUCCESS);
		//we get the size of the offset chuncks file	
		
		//
		/* Process all chunks */
		res = MPI_File_open(MPI_COMM_WORLD, file_tmp, MPI_MODE_RDONLY|MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, &fh_tmp1);
		assert(res == MPI_SUCCESS);
		
		res = MPI_File_open(MPI_COMM_WORLD, file_out, MPI_MODE_WRONLY|MPI_MODE_APPEND, MPI_INFO_NULL, &fh_out);
		assert(res == MPI_SUCCESS);
			
		//fprintf(stderr, "%s: rank %d : file_offset_sz= %zu \n", __func__, rank_num, file_offset_sz);
		MPI_Barrier(MPI_COMM_WORLD);
		if (file_r1 != NULL) {
			res = MPI_File_open(MPI_COMM_WORLD, file_r1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_r1);
			assert(res == MPI_SUCCESS);
		}		
		if (file_r2 != NULL) {
			res = MPI_File_open(MPI_COMM_WORLD, file_r2, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_r2);
			assert(res == MPI_SUCCESS);
		}		

		MPI_Barrier(MPI_COMM_WORLD);
		//fprintf(stderr, "%s: rank %d : start loop  \n", __func__, rank_num);
		buffer_r1 = buffer_r2 = NULL; seqs = NULL;
		coff = malloc(2 * sizeof(*coff));
		assert(coff != NULL);
		int chunck_num = 0;
		off_t offset_read = 0;
		off_t elmt_sz = 2*sizeof(size_t);
		off_t delta = proc_num * elmt_sz;		
		off_t previous_offset_chunck = 0;

		while (1) {
			char *p, *q, *e;
			int line_number;
			size_t reads_r1, reads_r2, reads;
			int64_t bases;

			/* Get current chunk offset and size ... */
			//res = MPI_File_read_shared(fh_tmp1, coff, 2, MPI_OFFSET, &status);
			// problem encountered with the read shared does not seem to work on NFS 
			// in replacement we do a a read_at
			offset_read = (rank_num * elmt_sz) + chunck_num * delta;
			if ( offset_read > file_offset_sz) break;
			res = MPI_File_read_at(fh_tmp1, offset_read, coff, 2, MPI_OFFSET, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_Get_count(&status, MPI_OFFSET, &count);
			assert(res == MPI_SUCCESS);
			//if ( coff[0] < previous_offset_chunck ) break;
			//if ( coff[1] == 0 ) break;
			//fprintf(stderr, "%s: rank %d : read offset[0] = %zu\n", __func__, rank_num, coff[0] );
			//fprintf(stderr, "%s: rank %d : read offset[1] = %zu\n", __func__, rank_num, coff[1] );
			if (count == 0) break;
			assert(count == 2);
			/* Read sequence datas ... */
			bef = MPI_Wtime();

			if (file_r1 != NULL) {
				buffer_r1 = realloc(buffer_r1, (size_t)coff[1]);
				assert(buffer_r1 != NULL);
				res = MPI_File_read_at(fh_r1, coff[0], buffer_r1, (int)coff[1], MPI_CHAR, &status);
				assert(res == MPI_SUCCESS);
				res = MPI_Get_count(&status, MPI_CHAR, &count);
				assert(res == MPI_SUCCESS);
				assert(count == (int)coff[1]);
				assert (*buffer_r1 == '@');
				
			}
			if (file_r2 != NULL) {
				buffer_r2 = realloc(buffer_r2, (size_t)coff[1]);
				assert(buffer_r2 != NULL);
				res = MPI_File_read_at(fh_r2, coff[0], buffer_r2, (int)coff[1], MPI_CHAR, &status);
				assert(res == MPI_SUCCESS);
				res = MPI_Get_count(&status, MPI_CHAR, &count);
				assert(res == MPI_SUCCESS);
				assert(count == (int)coff[1]);
				assert(*buffer_r2 == '@');
			}
			 chunck_num++;
			aft = MPI_Wtime();
			fprintf(stderr, "%s: read sequences (%.02f)\n", __func__, aft - bef);
			/* Count sequences ... */
			bef = MPI_Wtime();
			reads_r1 = reads_r2 = 0;
			if (file_r1 != NULL) {
				p = buffer_r1; e = buffer_r1 + coff[1];
				while (p < e) { reads_r1 += (*p++ == '\n') ? 1 : 0; }
			}
			if (file_r2 != NULL) {
				p = buffer_r1; e = buffer_r1 + coff[1];
				while (p < e) { reads_r2 += (*p++ == '\n') ? 1 : 0; }
			}
			reads_r1 /= 4; reads_r2 /= 4;
			assert(reads_r2 == 0 || reads_r1 == reads_r2);
			reads = reads_r1 + reads_r2; bases = 0;
			assert(reads <= INT_MAX);
			total_local_reads += reads;
			/* Parse sequences ... */
			seqs = malloc(reads * sizeof(*seqs));
			assert(seqs != NULL);
			if (file_r1 != NULL) {
				p = q = buffer_r1; e = buffer_r1 + coff[1]; line_number = 0;
				while (q < e) {
					if (*q != '\n') { q++; continue; }
					/* We have a full line ... process it */
					*q = '\0'; n = files * (line_number / 4);
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
						bases += seqs[n].l_seq;
						break;
					case 2: /* Line3: Ignored */
						assert(*p == '+');
						break;
					case 3: /* Line4: Quality */
						seqs[n].qual = p;
						break; }
					p = ++q; line_number++; }
			}
			if (file_r2 != NULL) {
				p = q = buffer_r2; e = buffer_r2 + coff[1]; line_number = 0;
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
						bases += seqs[n].l_seq;
						break;
					case 2: /* Line3: Ignored */
						assert(*p == '+');
						break;
					case 3: /* Line4: Quality */
						seqs[n].qual = p;
						break; }
					p = ++q; line_number++; }
			}
			aft = MPI_Wtime();
			xfprintf(stderr, "%s: parsed sequences (%.02f)\n", __func__, aft - bef);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] read %zu sequences (%ld bp)...\n", __func__, reads, (long)bases);

			/* Datas computation ... */
			bef = MPI_Wtime();
			mem_process_seqs(opt, indix.bwt, indix.bns, indix.pac, 0, (int)reads, seqs, pes0);
			aft = MPI_Wtime();
			fprintf(stderr, "%s: computed mappings (%.02f)\n", __func__, aft - bef);

			/* Write results ... */
			bef = MPI_Wtime();
			localsize = 0;
			for (n = 0; n < reads; n++) {
				/* Reuse .l_seq to store SAM line length to avoid multiple strlen() calls */
				seqs[n].l_seq = strlen(seqs[n].sam);
				localsize += seqs[n].l_seq; }
			assert(localsize <= INT_MAX);
			buffer_out = malloc(localsize);
			assert(buffer_out != NULL);
			p = buffer_out;
			for (n = 0; n < reads; n++) {
				memmove(p, seqs[n].sam, seqs[n].l_seq);
				p += seqs[n].l_seq;
				free(seqs[n].sam); }
			free(seqs);
			res = MPI_File_write_shared(fh_out, buffer_out, (int)localsize, MPI_CHAR, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_Get_count(&status, MPI_CHAR, &count);
			assert(res == MPI_SUCCESS);
			assert(count == (int)localsize);
			free(buffer_out);
			aft = MPI_Wtime();
			fprintf(stderr, "%s: wrote results (%.02f)\n", __func__, aft - bef);
			previous_offset_chunck = coff[0];

		}
		//fprintf(stderr, "%s: rank : %d : finish aligning = \n", __func__, rank_num);
		free(coff);
		free(buffer_r1);
		free(buffer_r2);

		res = MPI_File_close(&fh_tmp1);
		assert(res == MPI_SUCCESS);
		
		total_sz += tmp_sz;
		tmp_off +=  tmp_sz;	
		
		//we compute the total reads processed
		res = MPI_Reduce(&total_local_reads, &total_reads, 1,MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);		
		assert(res == MPI_SUCCESS);
		if (rank_num == 0)
			fprintf(stderr, "%s: total read processed  = %zu \n", __func__, total_reads);

	}//end if if (stat_r1.st_size == stat_r2.st_size)
	else {

		/*
		 * We are in the case the reads are trimmed
		 * Now each proc load the index and place it in shmem
		 * location
		 *
		 */

		if (proc_num < 2){

			fprintf(stderr, "For paired reads trimmed. The number of jobs must be greater than 1.\n"
					"Because one job cares for the forward reads and the other for the backward reads.\n"
					"This limitation should disappear in futur \n");

			res = MPI_Finalize();
			assert(res == MPI_SUCCESS);
			return 0;
		}

		fprintf(stderr, "in trimmed part\n");

		/***************************************************/

		/*
		 * We create 2 groups of jobs
		 *
		 * Say the even rank parse the backward
		 * and the odd parse the forward
		 *
		 * first part
		 *
		 * The odd ranks parse the forward file and keep track of the read number and the start offset of the read in the file
		 * The even ranks parse the backward file and keep track of the read number and the start offset of the read in the file
		 *
		 *
		 * second part
		 *
		 * Then even ranks trade the local number of read and local offset and build an index of the all reads
		 * with their offsets in the fastq file. Communication are done in the splitted communicator.
		 *
		 * third part
		 *
		 * we return in the global communicator to merge the 2 tables.
		 * The idea is to save in fd_tmp 3 columns
		 * read number | start offset read forward | start offset read backward
		 *
		 *
		 */


		size_t chunk_count = 0;
		int color = rank_num%2;

		MPI_Comm INTER_COMM;

		/*
		 * Split sequence files in chunks
		 */
		MPI_File mpi_fd_in1;
		bef = MPI_Wtime();
		char *file_to_read;
		if (color == 1){
			//odd ranks get the forward fastq
			file_to_read = file_r1;
		}
		else{
			//even ranks get the backward fastq
			file_to_read = file_r2;
		}
		MPI_File_open(MPI_COMM_WORLD, file_to_read,  MPI_MODE_RDONLY , MPI_INFO_NULL, &mpi_fd_in1);


		/*
		 * first we parse the buffer and see
		 * how many reads we have
		 */

		size_t read_number=0;
		assert(fd_in1 != -1);
		size_t *goff = NULL; //global offset contain the start offset in the fastq
		goff =(size_t*)calloc((size_t)(proc_num/2+1), sizeof(size_t));

		if (color == 1){

			//if rank is 1,3,5,7
			//We place file offset of each process to the begining of one read's line
			init_goff(goff, mpi_fd_in1, stat_r1.st_size, proc_num/2, rank_num);
			//we compute the new size according to the shift
			//We calculate the size to read for each process
			int ind = 0;
			ind = (int)floor(rank_num/2);
			size_t siz2read;
			siz2read = goff[ind+1]-goff[ind];
			size_t j=0;
			size_t off_in_file = goff[ind]; //Current offset in file sam
			buffer_r1 = (char*)calloc(siz2read+1,sizeof(char));
			buffer_r1[siz2read] = 0;
			MPI_File_read_at(mpi_fd_in1, (MPI_Offset)off_in_file, buffer_r1, siz2read, MPI_CHAR, MPI_STATUS_IGNORE);
			assert(buffer_r1 != MAP_FAILED);
			assert(*buffer_r1 == '@');
			p = buffer_r1 ;
			e = buffer_r1 + siz2read;
		}
		else{

			//if rank is even
			// 0,2,4,6,8,10
			//We place file offset of each process to the begining of one read's line
			init_goff(goff ,mpi_fd_in1, stat_r2.st_size, proc_num/2, rank_num);
			//we compute the new size according to the shift
			//We calculate the size to read for each process
			int ind = rank_num/2;
			size_t siz2read = goff[ind+1]-goff[ind];
			size_t j=0;
			size_t off_in_file = goff[ind]; //Current offset in file sam
			buffer_r1 = (char*)calloc(siz2read+1,sizeof(char));
			buffer_r1[siz2read] = 0;
			MPI_File_read_at(mpi_fd_in1, (MPI_Offset)off_in_file, buffer_r1, siz2read, MPI_CHAR, MPI_STATUS_IGNORE);
			assert(buffer_r1 != MAP_FAILED);
			assert(*buffer_r1 == '@');
			p = buffer_r1 ;
			e = buffer_r1 + siz2read;
		}

		//we count the number of lines
		int m = 0;
		char *t = p;
		char *u = p;
		int lines = 0;
		while (t < e){
			if (*t == '\n') lines++;
				t++;
		}
		size_t local_num_reads=lines/4;
		fprintf(stderr, "rank %d ::: number of reads = %zu \n", rank_num, local_num_reads);


		//we malloc a vector to hold offsets of the read
		//local_read_offset contains the ofset of the read
		//in the file
		size_t local_read_offsets[local_num_reads];
		lines = 0;
		if (color == 1){
			int ind = 0;
			ind = (int)floor(rank_num/2);
			local_read_offsets[0] =goff[ind];}
		else{local_read_offsets[0] =goff[rank_num/2];}

		int n = 1;
		size_t g=0;
		while (u < e && n < local_num_reads){
			if (*u == '\n') lines++;
			if (lines == 4) {
				local_read_offsets[n] = local_read_offsets[0] + g +1; //+1 for \n caracters
				n++; lines = 0;}
				u++;g++;
		}
		if (buffer_r1)
			free(buffer_r1);

		MPI_File_close(&mpi_fd_in1);
		//the roots in NEW_COMM compute the total number of reads
		size_t total_num_reads = 0;
		MPI_Reduce(&local_num_reads, &total_num_reads, 1,MPI_LONG_LONG_INT, MPI_SUM, 0,MPI_COMM_WORLD);
		//MPI_Reduce(&local_num_reads, &total_num_reads, 1,MPI_LONG_LONG_INT, MPI_SUM, 1,MPI_COMM_WORLD);
		total_num_reads = total_num_reads/2;
		MPI_Bcast(&total_num_reads, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

		//Now the roots of new communicator gather all the offsets vectors for R1 and R2
		//size_t *num_reads_by_proc = (size_t *)malloc(new_rank_num * sizeof(size_t));
		size_t num_reads_by_proc[proc_num];
		res = MPI_Allgather(&local_num_reads, 1, MPI_LONG_LONG_INT, num_reads_by_proc, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
		//we dispatch the previous result
		//

		if ( rank_num == 0 || rank_num ==1 )
			fprintf(stderr, "rank %d ::: total_num_reads = %zu \n", rank_num, total_num_reads);

		//malloc don't work problem in malloc_wrap (line 26)
		//size_t *all_reads_offsets = malloc(total_num_reads*sizeof(size_t));
		size_t all_reads_offsets[total_num_reads];
		all_reads_offsets[0]=0;
		assert(all_reads_offsets !=0);

		if (rank_num == 1){
					//first we copy local_read_num from 0
					// we compute displacement in both cases
					// odd and even ranks
					size_t disp[proc_num/2 + 1];
					assert(disp !=0);

					disp[0] = 0;

					for (m = 1; m < (proc_num/2 +1); m++)
						disp[m]= num_reads_by_proc[2*m-1];

					for ( m = 1; m < (proc_num/2 +1); m++)
						disp[m] = disp[m] + disp[m-1];
					//we copy for root rank
					for (m = 0; m < local_num_reads; m++){
						all_reads_offsets[m] = local_read_offsets[m];
					}

					for ( m = 1; m < proc_num/2; m++){

						size_t tmp[num_reads_by_proc[2*m +1]];
						MPI_Recv(tmp,
					             num_reads_by_proc[2*m+1],
					             MPI_LONG_LONG_INT,
					             2*m + 1,
					             0,
					             MPI_COMM_WORLD,
					             MPI_STATUS_IGNORE);

						size_t h=0;
						for (h = 0; h < num_reads_by_proc[2*m +1]; h++)
							all_reads_offsets[disp[m]+h] = tmp[h];
					}

					for ( m = 0; m < (total_num_reads - 1); m++){
						assert(all_reads_offsets[m] < all_reads_offsets[m + 1]);
					}
					assert(all_reads_offsets[total_num_reads-1] < stat_r1.st_size);

				}

		if (color == 1 && rank_num > 1){
			MPI_Send(local_read_offsets ,local_num_reads, MPI_LONG_LONG_INT, 1, 0, MPI_COMM_WORLD);
		}


		//each root in communicator NEW_COMM recieve all_reads_offsets
		if (rank_num == 0){

			// we compute displacement
			size_t disp[proc_num/2 + 1];
			assert(disp !=0);

			disp[0] = 0;
			for (m = 1; m < (proc_num/2 +1); m++)
				disp[m]= num_reads_by_proc[2*(m-1)];
			for ( m = 1; m < (proc_num/2 +1); m++)
				disp[m] = disp[m] + disp[m-1];

			//first we copy local_read_num from 0
			for (m = 0; m < local_num_reads; m++){
				all_reads_offsets[m] = local_read_offsets[m];
			}

			for ( m = 1; m < proc_num/2; m++){

				//size_t *tmp =malloc(sizeof(size_t)*num_reads_by_proc[2*m]);
				size_t tmp[num_reads_by_proc[2*m]];
				tmp[0] = 0;
				MPI_Recv(tmp,
		                num_reads_by_proc[2*m],
		                MPI_LONG_LONG_INT,
		                2*m,
		                0,
		                MPI_COMM_WORLD,
		                MPI_STATUS_IGNORE);

				size_t h=0;
				for (h = 0; h < num_reads_by_proc[2*m]; h++)
					all_reads_offsets[disp[m]+h] = tmp[h];

			}

			//FOR DEBUG

			for ( m = 0; m < (total_num_reads - 1); m++){
				assert(all_reads_offsets[m] < all_reads_offsets[m + 1]);
			}
			assert(all_reads_offsets[total_num_reads-1] < stat_r2.st_size);

			//fprintf(stderr, "rank %d ::: finish asserting \n", rank_num);
		}

		if (color == 0 && rank_num >0){
			MPI_Send(local_read_offsets ,local_num_reads, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);
		}

		//MPI_Barrier(MPI_COMM_WORLD);

		fprintf(stderr, "rank %d ::: finish copy offsets \n", rank_num);
		/*
		 *
		 * FOR DEBUG ::
		 *
		 * Verify the offset and print reads
		 * at some offsets

		if (rank_num == 0){


			size_t start=all_reads_offsets[1004];
			size_t end=all_reads_offsets[1005];

			char *buff = (char*)malloc((end-start)*sizeof(char));
			buff[end-start]=0;
			MPI_File_read_at(mpi_fd_in1, (MPI_Offset)start, buff, end-start , MPI_CHAR, MPI_STATUS_IGNORE);

			fprintf(stderr, "rank %d ::: Reads at offset %zu = %s \n", rank_num, start, buff);
			fprintf(stderr, "\n\n\n");
		}
		if (rank_num == 1){


			size_t start=all_reads_offsets[1004];
			size_t end=all_reads_offsets[1005];

			char *buff = (char*)malloc((end-start)*sizeof(char));
			buff[end-start]=0;
			MPI_File_read_at(mpi_fd_in1, (MPI_Offset)start, buff, end-start , MPI_CHAR, MPI_STATUS_IGNORE);

			fprintf(stderr, "rank %d ::: Reads at offset %zu = %s \n", rank_num, start, buff);
			fprintf(stderr, "\n\n\n");
		}
		*/

		/*
		 *
		 * In this part the root who cares for the forward (the rank 1)
		 * Compute chunk of reads. Each chunk has a size
		 * of approximately 100K pb.
		 * Then root send the index of offset to other root
		 *
		 */

		fprintf(stderr, "rank %d ::: compute chunk \n", rank_num );

		if (rank_num == 1){
			/*
			 * We compute number of chunks for forward reads
			 * and send it to
			 *
			 */
			size_t u=0;
			maxsiz = 100000;
			for (u = 0; u < total_num_reads; u++){
				//we for offsets multiple of maxsize
				if ( all_reads_offsets[u] > maxsiz){
					chunk_count +=1;
					maxsiz +=100000;
					//fprintf(stderr, "rank %d ::: chunk_count = %d \n", rank_num, chunk_count );
					//fprintf(stderr, "rank %d ::: position = %zu \n", rank_num, u );
					//fprintf(stderr, "rank %d ::: offset = %zu \n", rank_num, all_reads_offsets[u] );
				}
			}

			//we malloc a position vector
			//we will send to other root

			size_t position[chunk_count];
			position[0] = 0;
			chunk_count = 1;
			maxsiz = 100000;
			for (u = 0; u < total_num_reads; u++){
				//we for offsets multiple of maxsize
				if ( all_reads_offsets[u] > maxsiz){
					position[chunk_count] = u;
					chunk_count +=1;
					maxsiz +=100000;
				}
			}


			// now each root compute
			// coff an c.
			// coff contains the start offset of the chunks and the number of characters to
			// read

			coff = malloc(sizeof(MPI_Offset)*(chunk_count*2)); //hold the start offset of the chunk

			size_t last_chunk_start_offset = all_reads_offsets[position[chunk_count - 1]];
			size_t last_chunk_size = all_reads_offsets[total_num_reads -1] - all_reads_offsets[position[chunk_count - 1]];

			coff[0] = 0;
			coff[1] = all_reads_offsets[position[1]];
			size_t k=0;

			for (u = 1; u < (chunk_count - 1); u++){
				coff[2*u] = all_reads_offsets[position[u]];
				coff[2*u+1] = all_reads_offsets[position[u +1]] - all_reads_offsets[position[u]];
			}
			coff[2*chunk_count -2] = last_chunk_start_offset;
			coff[2*chunk_count -1] = last_chunk_size;


			//we send the number of chunks we found in forward files
			MPI_Send(&chunk_count, 1, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);

			//then we send the index of start offsets for each chunk
			MPI_Send(&position, chunk_count, MPI_LONG_LONG_INT, 0, 1, MPI_COMM_WORLD);

			fprintf(stderr, "rank %d ::: we send chunk count =%zu\n", rank_num, chunk_count);

		}

		if (rank_num == 0){
			//the root process recieve the number of position
			size_t u=0;
			MPI_Recv(&chunk_count,
					1,
					MPI_LONG_LONG_INT,
					1,
					0,
					MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);
			//size_t *position = malloc(sizeof(size_t)*chunk_count);
			size_t position[chunk_count];

			MPI_Recv(&position,
					 chunk_count,
					 MPI_LONG_LONG_INT,
					 1,
					 1,
					 MPI_COMM_WORLD,
					 MPI_STATUS_IGNORE);

			//now each root compute coff an c
			coff = malloc(sizeof(MPI_Offset)*(chunk_count*2)); //hold the start offset of the chunk

			size_t last_chunk_start_offset = all_reads_offsets[position[chunk_count - 1]];
			size_t last_chunk_size = all_reads_offsets[total_num_reads -1] - all_reads_offsets[position[chunk_count - 1]];

			coff[0] = 0;
			coff[1] = all_reads_offsets[position[1]];
			size_t k=0;

			for (u = 1; u < (chunk_count - 1); u++){
				coff[2*u] = all_reads_offsets[position[u]];
				coff[2*u+1] = all_reads_offsets[position[u +1]] - all_reads_offsets[position[u]];
			}

			coff[2*chunk_count -2] = last_chunk_start_offset;
			coff[2*chunk_count -1] = last_chunk_size;

	}

		/*
		 * FOR DEBUG
		 *
		 * Roots rank verify that we start at the right offsets
		 *

		if (rank_num == 0){

			int t=0;
			for (t = 0; t< chunk_count; t++){
			    size_t start=coff[2*t];
			    size_t end=start+5;

				char *buff = (char*)malloc((end-start)*sizeof(char));
				buff[end-start]=0;
				MPI_File_read_at(mpi_fd_in1, (MPI_Offset)start, buff, end-start , MPI_CHAR, MPI_STATUS_IGNORE);

				assert(buff[0] == '@');

				//fprintf(stderr, "rank %d ::: Reads at offset %zu = %s \n", rank_num, start, buff);
				//fprintf(stderr, "\n");
			}
		}
		if (rank_num == 1){
			int t=0;
			for (t = 0; t< chunk_count; t++){
			    size_t start=coff[2*t];
			    size_t end=start+5;

				char *buff = (char*)malloc((end-start)*sizeof(char));
				buff[end-start]=0;
				MPI_File_read_at(mpi_fd_in1, (MPI_Offset)start, buff, end-start , MPI_CHAR, MPI_STATUS_IGNORE);
				assert(buff[0] == '@');

				//fprintf(stderr, "rank %d ::: Reads at offset %zu = %s \n", rank_num, start, buff);
				//fprintf(stderr, "\n");
			}
		}
		*/


		//MPI_File_close(&mpi_fd_in1);

		//now we collect the vectors of offset for reads forward and backward
		//MPI_Barrier(NEW_COMM);
		MPI_Barrier(MPI_COMM_WORLD);
		/*
		 *
		 * Each Roots save offsets and size
		 * Save offsets+sizes in temp file
		 * PROBLEM:
		 *
		 * When file are going to read the 2 offsets there may be concurrency
		 * between access.
		 *
		 * A rank could read different offsets chunk.
		 *
		 *
		 */

		if (rank_num == 0){


			/*
			 * Rank 0 is going to write offsets and size for  forward and backward
			 */
			size_t g = 0;
			size_t *backward_offsets = malloc(chunk_count*2*sizeof(size_t));

			MPI_Recv(backward_offsets, chunk_count*2,
					MPI_LONG_LONG_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			// now we build a vector for all offset
			// size = chunk_count * 4
			size_t *all_offsets = malloc(chunk_count*4*sizeof(size_t));

			size_t k=0;
			for (g = 0; g < chunk_count*4; g++){

				//first we save forward offset
				all_offsets[g] = coff[2*k];
				all_offsets[g+1] = coff[2*k+1];
				all_offsets[g+2] = backward_offsets[2*k];
				all_offsets[g+3] = backward_offsets[2*k+1];
				g+=3;
				k++;
			}
			fprintf(stderr, "rank %d :::: Save offsets + sizes in temp file for even rank\n", rank_num);
			//we rename file tmp
			sprintf(file_tmp, "%s.offsets.tmp", file_out);

			res = MPI_File_delete(file_tmp, MPI_INFO_NULL);
			assert(res == MPI_SUCCESS || res == MPI_ERR_NO_SUCH_FILE || res == MPI_ERR_IO);
			res = MPI_File_open(MPI_COMM_SELF, file_tmp, MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &fh_tmp);
			assert(res == MPI_SUCCESS);
			res = MPI_File_write(fh_tmp, all_offsets, chunk_count*4, MPI_OFFSET, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_Get_count(&status, MPI_OFFSET, &count);
			assert(res == MPI_SUCCESS);
			assert(count == chunk_count*4);
			res = MPI_File_close(&fh_tmp);
			assert(res == MPI_SUCCESS);
			free(coff);
			MPI_File_close(&fh_tmp);

		}
		if (rank_num == 1){

			MPI_Send(coff,
					chunk_count*2,
					MPI_LONG_LONG_INT,
					0,
					0,
					MPI_COMM_WORLD);
			free(coff);

		}

		MPI_Barrier(MPI_COMM_WORLD);
		if (rank_num == 0)
			fprintf(stderr, "chunks and offsets saved\n");


		/*
		 * Now we dont' need split comm.
		 * we should MPI_Comm_free(&NEW_COMM);
		 * We return to MPI_COMM_WORLD
		 */

		/* Process all chunks */
		char file_tmp_off[PATH_MAX];
		sprintf(file_tmp_off, "%s.offsets.tmp", file_out);


		MPI_File fh_tmp_offsets; //for forward offsets

		/*
		 * Problem:
		 *
		 * A rank should access the same offsets chunk in 2 files.
		 *
		 */

		res = MPI_File_open(MPI_COMM_WORLD, file_tmp_off, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_tmp_offsets);
		assert(res == MPI_SUCCESS);

		res = MPI_File_open(MPI_COMM_WORLD, file_out, MPI_MODE_WRONLY|MPI_MODE_APPEND, MPI_INFO_NULL, &fh_out);
		assert(res == MPI_SUCCESS);

		if (file_r1 != NULL) {
			//fprintf(stderr, "rank %d ::: open file %s \n",rank_num, file_r1);
			res = MPI_File_open(MPI_COMM_WORLD, file_r1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_r1);
			assert(res == MPI_SUCCESS);
		}
		if (file_r2 != NULL) {
			//fprintf(stderr, "rank %d ::: open file %s \n",rank_num, file_r2);
			res = MPI_File_open(MPI_COMM_WORLD, file_r2, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh_r2);
			assert(res == MPI_SUCCESS);
		}

		buffer_r1 = buffer_r2 = NULL; seqs = NULL;
		MPI_Offset *coff_offsets = malloc(4 * sizeof(*coff));
		assert(coff_offsets != NULL);

		// here we loop until there's nothing to read
		// in the offset and size file
		while (1) {

			char *p, *q, *e;
			int line_number, line_number2;
			size_t reads_r1, reads_r2, reads;
			int64_t bases;

			/* Get current chunk offset and size for forward reads */
			res = MPI_File_read_shared(fh_tmp_offsets, coff_offsets, 4, MPI_OFFSET, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_Get_count(&status, MPI_OFFSET, &count);
			assert(res == MPI_SUCCESS);

			if (count == 0) break;

			assert(count == 4);
			/*
			fprintf(stderr, "rank :::: %d After open read shared coff_forward[0] = %zu ::: coff_forward[1] = %zu \n",
					rank_num, coff_offsets[0], coff_offsets[1]);

			fprintf(stderr, "rank :::: %d After open read shared coff_backward[2] = %zu ::: coff_backward[3] = %zu \n",
					rank_num, coff_offsets[2], coff_offsets[3]);
			*/


			/*
			 *
			 * Read sequence datas ...
			 *
			 */
			bef = MPI_Wtime();

			buffer_r1 = realloc(buffer_r1, (size_t)coff_offsets[1]+1);
			assert(buffer_r1 != NULL);

			buffer_r1[coff_offsets[1]]=0;
			res = MPI_File_read_at(fh_r2, coff_offsets[0], buffer_r1, coff_offsets[1], MPI_CHAR, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_Get_count(&status, MPI_CHAR, &count);
			assert(res == MPI_SUCCESS);
			assert(count == (int)coff_offsets[1] && *buffer_r1 == '@');


			buffer_r2 = realloc(buffer_r2, (size_t)coff_offsets[3]+1);
			assert(buffer_r2 != NULL);
			buffer_r2[coff_offsets[3]]=0;

			res = MPI_File_read_at(fh_r1, coff_offsets[2], buffer_r2, coff_offsets[3], MPI_CHAR, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_Get_count(&status, MPI_CHAR, &count);
			assert(res == MPI_SUCCESS);
			assert(count == (int) coff_offsets[3] && *buffer_r2 == '@');

			aft = MPI_Wtime();
			xfprintf(stderr, "%s: read sequences (%.02f)\n", __func__, aft - bef);

			/* Count sequences ... */
			bef = MPI_Wtime();
			reads_r1 = reads_r2 = 0;
			if (file_r1 != NULL) {
				p =buffer_r1; e = buffer_r1 + coff_offsets[1];
				while (p < e) { reads_r1 += (*p++ == '\n') ? 1 : 0; }
			}
			if (file_r2 != NULL) {
				p = buffer_r2; e = buffer_r2 + coff_offsets[3];
				while (p < e) { reads_r2 += (*p++ == '\n') ? 1 : 0; }
			}
			reads_r1 /= 4; reads_r2 /= 4;
			assert(reads_r1 == reads_r2);
			reads = reads_r1 + reads_r2; bases = 0;
			assert(reads <= INT_MAX);

			/* Parse sequences ... */
			seqs = malloc(reads * sizeof(*seqs));
			assert(seqs != NULL);
			if (file_r1 != NULL) {
				p = q = buffer_r1; e = buffer_r1 + coff_offsets[1]; line_number = line_number2 = 0;
				while (q < e) {
					if (*q != '\n') { q++; continue; }
					/* We have a full line ... process it */
					*q = '\0'; n = files * (line_number / 4);
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
						bases += seqs[n].l_seq;
						break;
					case 2: /* Line3: Ignored */
						assert(*p == '+');
						break;
					case 3: /* Line4: Quality */
						seqs[n].qual = p;
						break; }
					p = ++q; line_number++; }
				fprintf(stderr, "rank %d ::: Parse %d lines \n",rank_num, line_number );
			}
			if (file_r2 != NULL) {
				p = q = buffer_r2; e = buffer_r2 + coff_offsets[3]; line_number2 = 0;
				while (q < e) {
					if (*q != '\n') { q++; continue; }
					/* We have a full line ... process it */
					*q = '\0'; n = files * (line_number2 / 4) + 1;
					switch (line_number2 % 4) {
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
						bases += seqs[n].l_seq;
						break;
					case 2: /* Line3: Ignored */
						assert(*p == '+');
						break;
					case 3: /* Line4: Quality */
						seqs[n].qual = p;
						break; }
					p = ++q; line_number2++; }
				fprintf(stderr, "rank %d ::: Parse %d lines \n",rank_num, line_number );
			}

			assert(line_number == line_number2);

			//MPI_Barrier(MPI_COMM_WORLD);

			size_t g=0;
			for ( g =0; g < reads; g++){
				assert(seqs[g].name != NULL);
				assert(strlen(seqs[g].seq) ==  seqs[g].l_seq);
				assert(strlen(seqs[g].qual) ==  seqs[g].l_seq);
				/*
				 * FOR DEBUG
				 *
				if ((strlen(seqs[g].seq) !=  seqs[g].l_seq) || (strlen(seqs[g].qual) != seqs[g].l_seq)){

					fprintf(stderr, "rank %d ::: problem with read %zu \n",rank_num, g );
					fprintf(stderr, "rank %d ::: strlen(seqs[%zu].seq) = %zu \n",rank_num, g, strlen(seqs[g].seq) );
					fprintf(stderr, "rank %d ::: seqs[%zu].l_seq = %d \n",rank_num, g, seqs[g].l_seq );
					fprintf(stderr, "rank %d ::: strlen(seqs[%zu].qual) = %zu \n",rank_num, g, strlen(seqs[g].qual) );
				}
				*/
			}
			aft = MPI_Wtime();
			xfprintf(stderr, "%s: parsed sequences (%.02f)\n", __func__, aft - bef);
			if (bwa_verbose >= 3)
				fprintf(stderr, "rank %d ::: [M::%s] read %zu sequences (%ld bp)...\n",rank_num , __func__, reads, (long)bases);

			/* Datas computation ... */
			bef = MPI_Wtime();
			fprintf(stderr, "rank %d ::::  Call memseqs with count sequences %zu \n", rank_num, reads);

			mem_process_seqs(opt, indix.bwt, indix.bns, indix.pac, 0, (int)reads, seqs, pes0);

			aft = MPI_Wtime();
			fprintf(stderr, "%s: computed mappings (%.02f)\n", __func__, aft - bef);

			//MPI_Barrier(MPI_COMM_WORLD);
			/* Write results ... */
			bef = MPI_Wtime();
			localsize = 0;
			for (n = 0; n < reads; n++) {
				/* Reuse .l_seq to store SAM line length to avoid multiple strlen() calls */
				seqs[n].l_seq = strlen(seqs[n].sam);
				localsize += seqs[n].l_seq; }
			assert(localsize <= INT_MAX);
			buffer_out = malloc(localsize);
			assert(buffer_out != NULL);
			p = buffer_out;
			for (n = 0; n < reads; n++) {
				memmove(p, seqs[n].sam, seqs[n].l_seq);
				p += seqs[n].l_seq;
				free(seqs[n].sam); }
			free(seqs);
			res = MPI_File_write_shared(fh_out, buffer_out, (int)localsize, MPI_CHAR, &status);
			assert(res == MPI_SUCCESS);
			res = MPI_Get_count(&status, MPI_CHAR, &count);
			assert(res == MPI_SUCCESS);
			assert(count == (int)localsize);
			free(buffer_out);
			aft = MPI_Wtime();
			xfprintf(stderr, "%s: wrote results (%.02f)\n", __func__, aft - bef);


		} // end loop upon offsets


		if (buffer_r1) free(buffer_r1);
		if (buffer_r2) free(buffer_r2);
		if (coff_offsets) free(coff_offsets);
		if (goff) free(goff);

		res=MPI_File_close(&fh_tmp_offsets);
		assert(res == MPI_SUCCESS);
		//res = MPI_File_close(&fh_tmp); //closed before
		//assert(res == MPI_SUCCESS);

		fprintf(stderr, "rank %d :::: finish mappings for trimmed reads \n", rank_num);
	}

	free(opt);

	if (file_r2 != NULL) {
		res = MPI_File_close(&fh_r2);
		assert(res == MPI_SUCCESS);
	}
	if (file_r1 != NULL) {
		res = MPI_File_close(&fh_r1);
		assert(res == MPI_SUCCESS);
	}

	res = MPI_Win_free(&win_shr);
	assert(res == MPI_SUCCESS);

	res = MPI_File_close(&fh_out);
	assert(res == MPI_SUCCESS);

	bef = MPI_Wtime();
	res = MPI_Finalize();
	assert(res == MPI_SUCCESS);
	aft_total_run = MPI_Wtime();
	aft = MPI_Wtime();
	xfprintf(stderr, "%s: synched processes (%.02f)\n", __func__, aft - bef);
	fprintf(stderr, "%s: total run time (%.02f)\n", __func__, aft_total_run - bef_total_run);	
	return 0;
}
