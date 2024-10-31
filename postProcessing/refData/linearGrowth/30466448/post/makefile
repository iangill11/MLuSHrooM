#[
#[ Makefile for the MuSHrooM code.
#[ Darin Ernst, Manaure Francisquez.
#[
#[
#[ Users need to set the following environment variables:
#[
#[      FCOMPILER, FFTW_DIR, FFTW_LIB_openmp, ADIOS_DIR
#[
#[ This is done in setup.sh, which is invoked using ". setup.sh"  
#[ Below, defaults are used if the variables are not set.
#[


#[ Name of the program:
PROGNAME = ./shroom

#.Extract macros from header file.
useMPI := $(shell grep -e '\#define useMPI' shroomFLAGS.h | sed 's/\#define useMPI //' )
ioMode := $(shell grep -e '\#define ioMode' shroomFLAGS.h | sed 's/\#define ioMode //' )

all: $(PROGNAME).x

#[ Fortran compiler:
FC = $(FCOMPILER)

#[ Specify the location of FFTW and HDF5 here.
#[ These are generalized for the engage cluster (should work across accts)
FFTW_DIR ?= /usr
#[ These are for Jay Lang's Kali.
#[ FFTW_DIR ?= /home/linuxbrew/.linuxbrew/Cellar/fftw/3.3.8/
#[ These are for Mana's TheFarm.
#FFTW_DIR ?= /usr/local/Cellar/fftw/3.3.8/

#.User defined flag to use debug mode.
ifdef debug
   DEBUG_OPT = -g
   no_opt    =  1
else
   DEBUG_OPT =
endif

#.Flags used during testing and debugging.
#.THESE ARE FOR gfotran:
#CHECKS = #-Wall -Wrealloc-lhs -Wrealloc-lhs-all -g -ffpe-trap=invalid,zero,overflow -fcheck=bounds
#CHECKS = -Wall -ffpe-trap=invalid,zero,overflow -mcmodel=medium
#CHECKS = -Wall -ffpe-trap=invalid,zero,overflow -fcheck=bounds
#.THESE ARE FOR ifort:
CHECKS = #-fpe0 -no-prec-div -check bounds -warn all
#.Optimization flags. 
OPTIM  = -O3

ifdef openmp
   OPENMP = -openmp
else
   OPENMP =
endif

#.FFTW library flags, paths and library files.
ifeq ($(useMPI),0)
   FFTWlib = -lfftw3 
else
   #.For running parallel jobs use the following:
   FFTWlib = -lfftw3 -lfftw3_omp -lfftw3_mpi -lm
endif
FFTW = -L $(FFTW_DIR)/lib -I $(FFTW_DIR)/include $(FFTWlib) -Wl,-rpath -Wl,$(FFTW_DIR)/lib

ifeq ($(ioMode),1)
   #.ADIOS library flags, paths and library files.
   ADIOS = ${ADIOS_INC} ${ADIOS_FLIB} -L $(ADIOS_DIR)/lib -I $(ADIOS_DIR)/include
else
   ADIOS = 
endif

ifeq ($(ioMode),1)
   #.MPI library flags, paths and library files.
   MPI = -L ${MPI_DIR}/lib -I ${MPI_DIR}/include
else
   MPI = 
endif

#.List of modules in shroom.F90. Delete them after building executable.
mods = parameters fields MPItools IOtools FFTmodule adaptLive initialize terminate nonlinearities timeSteppers


FCflags  = $(OPTIM) $(CHECKS) $(DEBUG_OPT) -lpthread -fopenmp
LIBnFLAG = $(FCflags) $(MPI) $(FFTW) $(OPENMP) $(ADIOS) 

clean: $(PROGNAME).x
	rm -f $(addsuffix .mod, $(mods))

$(PROGNAME).x: $(PROGNAME).F90
	$(FC) -o $(PROGNAME).x $(PROGNAME).F90 $(LIBnFLAG)

$(PROGNAME):
	echo 'Error: You need the .x extension'



