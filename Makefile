# Camille Couprie 
# october 2009

OBJDIR	= objects
PINKDIR	= PINK
CSDIR	= CSparse
ARGV = argv
PWSRC = src
VPATH = ${PINKDIR}:${CSDIR}:${ARGV}

CSINCLUDE = -I${CSDIR}/Include 
PINKINCLUDE = -I${PINKDIR} 
ARGVINCLUDE = -I${ARGV} 
PWINCLUDE = -Iinclude
INCL = ${CSINCLUDE} ${PINKINCLUDE} ${ARGVINCLUDE} ${PWINCLUDE}

OBJ=	objects/larith.o \
	objects/ccsort.o \
	objects/cccodimage.o \
	objects/gageodesic.o \
	objects/mccodimage.o \
	objects/mcimage.o \
	objects/mcindic.o \
	objects/mclifo.o \
	objects/random_walker.o \
	objects/lMSF.o \
	objects/MSF_RW.o \
	objects/mcrbt.o \
	objects/union_find.o \
	objects/image_toolbox.o \
	objects/cs_lu.o \
	objects/cs_lusol.o \
	objects/cs_malloc.o \
	objects/cs_util.o \
	objects/cs_multiply.o \
	objects/cs_compress.o \
	objects/cs_lsolve.o \
	objects/cs_scatter.o \
	objects/cs_cumsum.o \
	objects/cs_sqr.o \
	objects/cs_ipvec.o \
	objects/cs_amd.o \
	objects/cs_permute.o \
	objects/cs_transpose.o \
	objects/cs_counts.o \
	objects/cs_add.o \
	objects/cs_etree.o \
	objects/cs_leaf.o \
	objects/cs_fkeep.o \
	objects/cs_tdfs.o \
	objects/cs_usolve.o \
	objects/cs_spsolve.o \
	objects/cs_post.o \
	objects/cs_reach.o \
	objects/cs_dfs.o \
	objects/argv.o  
#	objects/cs_print.o \
#	objects/cs_norm.o \
#	
FLAGS = -Wall # -pg 
CXX = g++ 
CC = gcc

OPTIMISE=-O4
WARNINGS=-Wall #-Werror
CXXFLAGS= ${WARNINGS} ${OPTIMISE} -Wno-deprecated ${FLAGS}
CFLAGS = ${WARNINGS} ${OPTIMISE} ${FLAGS}

all:	${CS} ${OBJ}
	${MAKE} powerwatsegm.exe 

# make with the Intel compiler, a very good complement to gcc/g++
# much more efficient and with better diagnoses
intel: 
	${MAKE} CC=icc CXX=icpc all

debug:
	${MAKE} OPTIMISE='' all

# remove the asserts and the debug information
production:
	${MAKE} DEBUG='' FLAGS="-DNDEBUG" all

clean:
	rm -f *.exe; rm -f *~; rm -f $(OBJ);  #rm -f overlay*; rm -f mask*;  

powerwatsegm.exe: ${PWSRC}/powerwatsegm.c  $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCL) ${PWSRC}/powerwatsegm.c $(OBJ) -o powerwatsegm.exe 

$(OBJDIR)/%.o:	${PWSRC}/%.c 
	$(CXX) $(CXXFLAGS) ${INCL} -c $< -o $@ 

$(OBJDIR)/%.o:	${CSDIR}/Source/%.c ${CSDIR}/Include/cs.h
	$(CC) $(CFLAGS) ${CSINCLUDE} -c $< -o $@

$(OBJDIR)/%.o:	${PINKDIR}/%.c
	$(CXX) $(CXXFLAGS) ${PINKINCLUDE} -c $< -o $@

$(OBJDIR)/argv.o:	${ARGV}/argv.c ${ARGV}/argv.h ${ARGV}/argv_loc.h
	$(CC) $(CCFLAGS) ${ARGVINCLUDE} -c $< -o $@

$(OBJDIR)/lMSF.o:	${PWSRC}/lMSF.c ${PINKDIR}/mcimage.h ${PINKDIR}/mccodimage.h  
	$(CXX) $(CFLAGS) ${PINKINCLUDE} ${PWINCLUDE} -c ${PWSRC}/lMSF.c -o $@


