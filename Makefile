#
#  Makefile for CRUMP 
#

# Compiler and flags
FC = gfortran
FFLAGS = -O3 -fopenmp -std=f2003
FSRC = crump_mod.f90 

# Object files
OBJS = $(FSRC:.f90=.o) 

# Name of library
LIBNAME=crumpf.so

$(LIBNAME)	: $(OBJS)
	$(FC) $(FSRC) $(FFLAGS)  -fPIC -shared -o $(LIBNAME)

%.o : %.f90
	$(FC) $< -o $@ -c $(FFLAGS)      

clean :
	-\rm $(LIBNAME) *.o  *.mod 

