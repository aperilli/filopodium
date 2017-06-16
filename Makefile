SEC = second.o

#    GNU
#F77  = /usr/local/pgi/linux86/bin/pgf77 -O3 -r8
#F77 = gfortran -O3 -w -g
#F77 = gfortran -O3 -w #-ftree-vectorize
#F77 = gfortran -g  -C
#LIB = -L/sw/lib/ -llapack

####   XLF   
#F77 = /opt/ibmcmp/xlf/8.1/bin/f77 -O2 -WF,-DOSX -qmaxmem=-1 -qdpc -qextname 
#F77 = /opt/ibmcmp/xlf/8.1/bin/f77 -g -C -qdpc -qmaxmem=-1 #-qextname 
#F77 = /opt/ibmcmp/xlf/8.1/bin/f77 -O3 -WF,-DOSX -qmaxmem=-1 -qdpc #-qextname 
#XTRALIB = -Wl,-framework -Wl,Accelerate 

#Linux INTEL
#F77 = ifort -p -g -check bounds -w -r8
#F77 = ifort -g  -w -check bounds -r8
#F77 = mpif90 -g  -w -check bounds -fpe0 -r8
#F77 = mpif90 -g  -w -C 
#F77 = ifort -g -C -w -r8
#F77 = ifort -O3 -w -r8
#F77 = mpif90 -w -O3

F77 = mpiifort 
CC = mpiicc
F77_1 = ifort -fpe0 -r8 

#F77 = ${HOME}/bin/mympif77 -WF,-DOSX -qmaxmem=-1
#CC = cc -g -w -r8
#CC = cc -O3 -w 


COMPOPTS = -r8 -w -fpp #-Vaxlib -DPOINTER_SIZE=8 #-mcmodel medium -shared-intel   # INTEL
FFLAGS = -w -O3 -r8
#FFLAGS = -w -g -C -r8
#XTRALIB = 
RNGLIB = -L${HOME}/lib/sprng/lib -llcg64
INC = ${HOME}/lib/sprng/include
XTRALIB = #-L/usr/mpi/intel/openmpi-1.4.3/lib64 -lmpi

#COMPOPTS = -O3 -qarch=G5 -qtune=G5 -qhot -qipa -qstrict
COMPOPTS := $(COMPOPTS) -O3 -xMIC-AVX512 -fma -align array64byte -finline-functions -dynamic  #-mcmodel=large -shared-intel #-fpe0
COMPOPTS_C := $(COMPOPTS_C) -O3 -xMIC-AVX512 -fma -finline-functions -dynamic  # -mcmodel=large -shared-intel

DEPEND = Makefile

OBJS = st_mpi.o common_variables.o parallel_mpi.o compute_energy.o products.o filopod_wrapper.o \
       filopodium_mpi.o sprng_mpi.o ipickoff.o intread.o \
       read_input.o init_filaments.o compute_ekin.o compute_forces.o \
       randomgaussian.o integrate_velocity.o integrate_position.o integrate_filament.o \
       attempt_reaction.o polymerization.o depolymerization.o \
       integrate_wall.o integrate_Lwall.o integrate_vel_wall.o

.SUFFIXES:.f90 .c .inc .cm .o .mod

#.f.o:
#	$(F77) -c $<

ALL = rng_p fil_paral

fil_mpi : $(ALL) 

rng_p : 
	$(CC) $(COMPOPTS_C) -I$(INC) -c sprng_mpi.c

fil_paral: $(DEPEND) $(OBJS)
	$(F77) $(COMPOPTS) $(XTRALIB) -o fil_mpi $(OBJS) $(RNGLIB) 

%.o : %.f90
	$(F77) -c $<

clean:
	rm -f *.o fil_mpi

realclean:
	(cd LIB; make clean)
	(cd ANALYZE; make clean)
	(cd ANALYZE/KVECT; make clean)
	rm -f *.o fil_mpi

$(OBJS) : $(DEPEND)

