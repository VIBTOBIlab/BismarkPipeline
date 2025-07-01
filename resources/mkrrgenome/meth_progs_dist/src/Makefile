# make file for cpgis
#
# Configure the following to suit your setup..
#
# edit these to reflect actual library and include locations..
HDRLOC = ../include/

#
# define compile and load commands here.
CFLGS =
CC = cc
CCLDFLGS =
CINCLUDES = -I../include
CCLIBS = -lm
CCMP = $(CC) -c $(CFLGS) -o $@ $(CINCLUDES) $*.c
CCLD = $(CC) $(CCLDFLGS) -o $@ $@.o $(CCLIBS) $(SYSLIBS)
OBJS = bas_fns.o sqfl_fns.o sqmat_fns.o wlu_fns.o rbs_fns.o mrg_fns.o
C_OBJS = cpgis.o $(OBJS)
FSMOBJS = $(OBJS) fsm_ops.o

# dependency rules

EXES = rmapbsbed2cpg bin_cnts mkrrgenome scan_cpg_depth cleanadaptors

all: $(EXES)

rmapbsbed2cpg: rmapbsbed2cpg.o $(OBJS)
	       $(CCLD) $(OBJS)

rmapbsbed2cpg.o: rmapbsbed2cpg.c rmapbsbed2cpg.h mrg_fns.h rbs_fns.h
		 $(CCMP)	       

bin_cnts: bin_cnts.o $(FSMOBJS)
	  $(CCLD) $(FSMOBJS)

bin_cnts.o: bin_cnts.c bin_cnts.h
	    $(CCMP)

mkrrgenome: mkrrgenome.o $(FSMOBJS)
	    $(CCLD) $(FSMOBJS)

mkrrgenome.o: mkrrgenome.c bin_cnts.h rmapbsbed2cpg.h
	      $(CCMP)

scan_cpg_depth: scan_cpg_depth.o $(FSMOBJS)
		$(CCLD) $(FSMOBJS)

scan_cpg_depth.o: scan_cpg_depth.c
		  $(CCMP)

cleanadaptors: cleanadaptors.o $(FSMOBJS)
	  $(CCLD) $(FSMOBJS)

cleanadaptors.o: cleanadaptors.c
	    $(CCMP)

bas_fns.o: bas_fns.c $(HDRLOC)bas_fns.h
	$(CCMP)

sqfl_fns.o: sqfl_fns.c $(HDRLOC)sqfl_fns.h
	$(CCMP)

wlu_fns.o: wlu_fns.c $(HDRLOC)wlu_fns.h
	$(CCMP)

sqmat_fns.o: sqmat_fns.c $(HDRLOC)sqmat_fns.h
	$(CCMP)

fsm_ops.o: fsm_ops.c $(HDRLOC)fsm_ops.h
	$(CCMP)

rbs_fns.o: rbs_fns.c rbs_fns.h
	$(CCMP)

mrg_fns.o: mrg_fns.c mrg_fns.h
	$(CCMP)

clean:
	-rm *.o

CLEAN: clean
	-rm $(EXES)

debug:
	$(MAKE) all CFLGS=-g

mallocdebug: clean
	$(MAKE) all CFLGS="-g -DMALLOCDBG"
