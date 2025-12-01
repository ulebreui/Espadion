#Name of the executable
EXEC=espadion

#Compiler
COMP=GNU

#Debug option
DEBUG=0


#############################################################################
# Fortran compiler options and directives
#############################################################################
# GFORTRAN
ifeq ($(COMP),GNU)
F90 = gfortran
ifeq ($(OPENMP),1)
F90 = gfortran -fopenmp
endif
   FFLAGS +=  -x f95-cpp-input -O3 -frecord-marker=4 -ffree-line-length-none $(DEFINES)
ifeq ($(DEBUG),1)
   FFLAGS += -g -fbounds-check -fbacktrace -Wuninitialized -Wall
   FFLAGS += -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan
   #FFLAGS += -finit-real=nan
endif
endif
#INTEL
ifeq ($(COMP),INTEL)
F90 = ifort
ifeq ($(OPENMP),1)
F90 = ifort -qopenmp
endif
FFLAGS = -cpp $(DEFINES) -DWITHOUTMPI
FFLAGS += -O3 -fp-model source
endif

#############################################################################m
# Sources directories are searched in this exact order
BINDIR = .
VPATH  = ../:../f90
#############################################################################
# All objects
MODOBJ = precision.o ionisation_commons.o
MODOBJ+= io.o
#############################################################################
all:
	for file in $(MODOBJ); do unset MAKELEVEL ; $(MAKE) $$file; done
	$(MAKE) espadion
espadion: $(MODOBJ) espadion.o
	   $(F90) $(MODOBJ) espadion.o -o ./$(EXEC)

%.o:%.F
	$(F90) $(FFLAGS) -c $^ -o $@ $(LIBS_OBJ)
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@ $(LIBS_OBJ)
%.o:%.f95
	$(F90) $(FFLAGS) -c $^ -o $@ $(LIBS_OBJ)
FORCE:
#############################################################################
clean :
	rm -f *.o *.$(MOD) *__genmod.f90 *mod espadion
#############################################################################
