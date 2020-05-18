The project has moved here:

https://github.com/bioinfo-pf-curie/mpiBWA

Your are in the master branch of the mpiBWA project.

We recommand you to go to the Experimental branch (git checkout Experimental + git pull). 
The experimental branch has full reproducibility and a better scaling.
This branch will become the master in futur.

Other branches (master, FULLMPI, LAZYCHUNK ) are history.

Release notes
------------

Release 1.0 from the 10/07/2018<br />

Changes in Experimental branch 

1) Fix a bug during the mapping in shared memory of the reference genome
This bug didn't appear with openMPI version but Intel compiler complains.

2) Creation of a google group <br />
	https://groups.google.com/forum/#!forum/hpc-bioinformatics-pipeline <br />

3) To improve performances on Lustre file system removed the “suid” mount option (rw,nosuid,flock,lazystatfs) <br />
With Beegfs set the "flock" to "on" for reproducibility.<br />

Release 1.0 from the 30/04/2018 <br />

Changes in Experimental branch <br />

1) Add support for trimmed reads.  <br />
2) Be aware of the flock mode on parallel file system (Lustre, beegfs): flock must be on.
Otherwise reproducibility is not guaranteed. <br />

Release 1.0 from the 21/03/2018 <br />

1) Fix an inversion in file handle (line 1039 and 1056) <br />

Release 1.0 from 23/12/2017

changes in Experimental branch <br />

1) 100% reproducibility with the pipeline control (bwa mem -t 8) <br />
2) remove memory limits now offset are computed on 1gb buffer <br />
3) remove some memory leaks <br />
4) optimization of code <br />

Release 1.0 from 04/12/2017

Changes in the branch Experimental.

1) Fix a memory leak.

Release 1.0 from 01/12/2017

Changes in the branch Experimental. <br />
Fix the mpi_read_at buffer limit size.<br />


Release 1.0 from 29/11/2017

Add a new branch called Experimental. <br />
Warning: This is experimental work do not use in production. But test it and send us reports.<br />

Rationnal:<br />

The master and FULLMPI branches are made for full reproducibility (independant to the number of jobs) and accuracy but they reach the Amdah'ls law point. 
Indeed the locking file RMA implementation (the serialization when computing offsets) introduces a bottle neck we are not able to overpass with RMA technics. <br />
 
This is why we have implemented a new algorithm. Now a lot more master jobs are responsible for chuncking the data the way bwa does and present it to bwa-mem aligner. 
And instead of doing it linearly on the fastq now they do it independently and in parallel. With a little inter communication they adjust the chunks offets and sizes. 
This method removes the serialization bottle neck. <br />

As in the master and FULLMPI branches we obtain a full reproducibility and with a better efficiency and scalability. <br />

When testing this branch make sure the total number of jobs you take is (master jobs) * 8. <br />
8 is the number of aligner threads used by bwa-mem. <br />
According to bwa-mem policy all chunks are 10e6 by the number of threads nucleotide bases big. <br />

Remark: The initial buffer of each master jobs is limited to 2gb (due to mpi_read_at buffer size). <br />
This version does not work on trimmed reads. <br />

First results test on broadwell. <br />

Sample: <br />

NA12878 Illumina 300X 2x150 WGS from GIAB chinese trio.

The alignement is done with 352*8 = 2816 cpu<br />
352 = (Forward Fastq size in gb) / 2g,  8 is thenumber of bwa-mem aligner jobs (8 per master jobs)<br />

MPI parameters :<br />
mpi_run -n 352 -c 8<br />
or : <br />
MSUB -n 352<br />
MSUB -c 8<br />

alignment time: 26 mn<br />
Time to compute chunks: 8s<br />

Reproducibility: the pipeline has been tested tested with 5632 and 2816 cpu results are the same. <br />

Next step:<br />

1) remove initial buffer size limitation<br />
2) support trimmed reads<br />
3) need tests on low throughput infrastructures<br />

Release 1.0 from 21/11/2017

Changes in  LAZYCHUNCK branch:

1) forget a end condition

Release 1.0 from 20/11/2017

Changes in  LAZYCHUNCK branch:
1) remove memory leaks
2) update results section
3) test with 600Gb (x2) fastq files and 100 jobs: ok. (20 mn to approximate chuncks offset) 
4) reproducibility with constant number of jobs: ok.
 
Notes:
If want 100% reproducibility whatever the number of jobs use the master branch or the FULLMPI branch.
But if you want go faster use LAZYCHUNK branch. 
 
Release 1.0 from 17/11/2017

Changes in  LAZYCHUNCK branch

1) Improvement of the algorithm for window sizing.
2) Change the number of bases for the estimation line 541. 

Release 1.0 from 15/11/2017

Changes in LAZYCHUNK branch

1) Due to a mmap in the previous version the virtual memory may be big when the fastq file is large.
Some scheduler like Torque doesn't like this. We review the algorithm of computing chuncks.
first we divide the file in window of size (number of jobs) * 1Gb. On each window jobs approximate chuncks of 10 mega bases.
this way the virtual memory stay low. 
 
2) We also saw problem on some architecture with MPI file read shared we replace it with MPI file read at. 

3) We have also other projects of parallelization (sorting, marking duplicates, clustering... ) and we are looking for people willing to help us for developments and tests. Don't hesitate to contact me (frederic.jarlier@curie.fr). 


Release 1.0 from 30/06/2017

Major changes:
=======
28/11/2017
>>>>>>> 884482b956169fabe0e8696021fc0092c403d86f

First release

We have implemented a new algorithm.
In this algorithm you have master jobs and aligner jobs.
Master jobs are responsible for computing offset and chunk sizes.
This way all the chuncks have the same number of bases. 
Then master jobs send chunks to bwa-mem and are aligned.
Finally master jobs write ni the result sam file.

the total number of jobs muste be: (number of master) * 8 
8 is the number of threads used by bwa-mem.

First result test on broadwell of TGCC

Condition: Due to mpi_read_at buffer limit of 2g. 
You should limit the initial buffer read to 2gb per master jobs.
Futur release will remove this condition. 

Tested on the SRR2052 WGS with 352*8 cpu 
alignment time: 26 mn
Time to compute chunks: 8s
Scalability : ok
Reproducibility : ok 
Command line :mpi_run -n 352 -c 8
or : MSUB -n 352
MSUB -c 8

This version does not support trimmed read yet.

Installation
---------

git clone https://github.com/fredjarlier/mpiBWA.git 
git checkout Experimental
git pull
.export PATH=/PATH_TO/automake-1.15/bin:/PATH_TO/autoconf-2.69/bin:$PATH
 ./configure CC=/PATH_TO/mpicc
make && make intall

Requirements
------------

You need a C compiler as required for classic BWA program.

You need to install a version of openMPI. (see: https://www.open-mpi.org/)

You need a mpi compiler too. to check your mpi installation tell in a command window whereis mpirun, normally it is installed in /usr/bin/mpirun.

You need a low latency network and a parallel file system to use this branch

Your reads should be paired.

Options
------

All the BWA-MEM options are available in this version, of course according to the bwa release wrapped. 

Known issues:
-------------

Primary hits are reproduced between the serial version and the parallel but you can see differences in mapping position for alternate contigs. 
This problem stems from the randomization of multi-hits reads. When running with the same number of MPI jobs alternative positions are reproduced but when the number of jobs varies the positions can switch for secondary alignments.

How to integrate further version
--------------------------------

This version of mpiBWA has been build with 0.7.15 BWA version.
To integrate the 0.7.+:

1) git clone the 0.7.+ of BWA.

2) in the folder of bwa copy-pass the following function from mpiBWA:

makefile
main_parallel_version.c 
pidx.c

then make.

Compilation 
-----------

You need automake 1.15 for the installation.
You can install automake and autoconf in differents directories and export the path like this:
export PATH=../automake-1.15/bin:../autoconf-2.69/bin:$PATH

Download from git. In the folder mpiBWA type:
./configure && make install && make

or for distribution:
make dist
tar xzf .tar.gz
cd pbwa7-1.0
./configure && make install && make

for passing mpi path:
./configure CC=mpi_bin_path
add --prefix in configure if you need 

Results are 2 executables pidx and pbwa7.

Build a reference
-----------------
After the creation of the reference file with BWA, you need to build a mapped reference genome. 
To do that: pidx my_ref.fa (where my_ref.fa has been build with BWA).
Pidx build a reference with the extension .map (my_ref.fa.map). 
This reference will be mapped in share memory.

The pidx needs the following files my_ref.fa.sa, my_ref.fa.bwt, my_ref.fa.ann, my_ref.fa.pac, my_ref.fa.amb to construct the fa.map.
 Those files are generated with bwa index. It works also with the fasta.

How to manage the multithreading
----------------

You have the option to use the multithreading options. 
To do that you fix the number of nodes with mpirun and then the number of threads per nodes with the -t.

Examples:

1 server with 8 threads per server:
mpirun -n 1  pbwa7 mem -t 8... 

2 servers with 10 threads per servers
mpirun -n 2  pbwa7 mem -t 10... 

And so on.

Example of command lines
---------------------

Example of bash to run the pbwa7 with openmpi 1.10.7

We launch 30 master jobs with 8 threads each it make a total of 240 jobs 

ERROR=$PATH/error.txt
SAM=$PATH/mysample.sam
TASKS=30
PPN=8
mpirun -np $TASKS --map-by ppr:$PPN:socket:pe=$PE $PBWA mem -t 8 -o $SAM $REF $FASTQ1 $FASTQ2 &> $ERROR

Results
-------

Here results wo obtain from test realized with TGCC (Très Grand Centre de Calcul - Bruyères le Chatel - France).

![img](Results_TGCC_Broadwell.jpg)

Remarks
-------

1) Do not type the .map extension when you give the reference to pbwa7.

2) Hybrid mode is possible with -t options. MPI rank fix the number of servers and -t options the number of threads per job. 

3) For Lustre or parallel file system users: If you intend to run the mpiSort after the alignment you can tell the striping of the results. 
This is done according to the striping you set in the mpiSort program with the Lustre "lfs setstripe" command (lfs setstripe -c -1 -s 2gb .). 
The "-c -1" option tells Lustre to use all the file system servers and and "-s 2gb" is the size of contigues data blocks. 
For speed purpose reading commands (particularly MPI commands) are aligned on those blocks.

Notes: The striping (like with Hadoop) is the way your data is distributed among servers. This technic accelerates access files and meta file information.

4) Make sure you have the same version of MPI for compiling and running.

5) From our experiences some scheduler do not interpret shared memory and jobs ends up stuck for overloading memory.
In this case take less than 1gb per job per threads (if you take 8 threads). This is approximately the total reference plus the chunk size.

Questions:
--------

We have created a google forum for questions.

https://groups.google.com/forum/#!forum/hpc-bioinformatics-pipeline

Feel free to ask any question

Future work:
----------

1) Overlap writing, reading and aligning.

References:
---------

This work is based on the original bwa aligner written by Li et al.

Li H. and Durbin R. (2010) Fast and accurate long-read alignment with Burrows-Wheeler transform. Bioinformatics, 26, 589-595. [PMID: 20080505] 

Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v1 [q-bio.GN]

Latham R. et al. (2007)  Implementing MPI-IO Atomic Mode and Shared File Pointers Using MPI One-Sided Communication Authors

-------

The program has been developed by 

Frederic Jarlier from Institut Curie and 
Nicolas Joly from Institut Pasteur

and supervised by

Philippe Hupe from Institut Curie
