# Makefile

FC = gfortran
FFLAGS = -O2 -fbacktrace -fbounds-check

.SUFFIXES :
.SUFFIXES : .f90 .o

default: isotxs_io.exe

.f90.o :
	$(FC) $(FFLAGS) -c $*.f90

isotxs_io.exe: isotxs_io.o
	$(FC) $(FFLAGS) isotxs_io.o -o $@

clean:
	rm -f core *.o *.mod *.exe

query:
	cvs -n update -P

update:
	cvs update -P