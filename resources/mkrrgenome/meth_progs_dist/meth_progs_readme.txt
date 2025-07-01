	  Programs and Scripts for Bisulphite Sequence Data

Peter A. Stockwell
Dept of Biochemistry,
University of Otago,
Dunedin, New Zealand.

peter.stockwell@otago.ac.nz

----------------------------------------------------------------------


rmapbsbed2cpg
bin_cnts
scan_cpg_depth
mkrrgenome
cleanadaptors

illum2fasta.awk
mk4to1lines.awk
bismmethex2list.awk
tidyrrnams.awk

These programs have been written in the course of work on differential
CpG methylation of the human genome and represent a series of tools
for preparing, modifying and analysing these data.  The work has
particularly focussed on reduced representation (RR) bisulphite
sequencing (RRBS) as in Meissner, et al., 2008, Nature, 454, 766-770,
Smith, et al., 2009) Methods, 48, 226-232 and Gu, et al., 2011, Nature
Protocols, 6, 468-481.  This software should, none-the-less, have
wider application than that.

This is research software so that it is not necessarily easy to use,
although it should work correctly as intended.  It works from a Unix
(MacOS X) or Linux command line interface: the notes below describe the
functionality and format of intermediate files where appropriate.
This code has been generated and tested on a MacOS X system (10.6)
using gcc v4.2.1 and on various Linux platforms (RedHat, Centos,
Fedora, Ubuntu) but is written to compile and run on any appropriate
C compiler and environment.  The size of files and data required for
this work will generally require a 64 bit environment.  Awk scripts
have been developed for the Gnu AWK (gawk) distributed with MacOS X
but, again, they should run in other comparable environments.

----------------------------------------------------------------------
		 Unpacking and building instructions:

Distribution: source code is encorporated into a self-extracting
shell archive which can be run with /bin/sh, or equivalent, to
generate the complete sources and Makefile.  The distributed form is a
file meth_progs_dist.shar which is unpacked by:

/bin/sh meth_progs_dist.shar

generating a directory meth_progs_dist containing two directories:
include, containing some of the generic header files and src
containing the rest of the source code and Makefile.  The set can be
built with the following commands:

cd meth_progs_dist/src
make

or

make CFLGS=-O3

to optimise the code.  Some warnings may appear but these can be
ignored.  The executables will be in the src directory - no install
targets are provided, the completed executables (rmapbsbed2cpg
bin_cnts scan_cpg_depth mkrrgenome cleanadaptors) should be manually
copied to an appropriate directory (e.g. /usr/local/bin).  The awk
scripts are in the src directory and may similarly be copied somewhere
useful.

Documentation for the set is present in the file progs_doc.txt which
resides in the top level of the directories.
