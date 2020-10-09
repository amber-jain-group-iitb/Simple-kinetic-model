MKLROOT=/opt/intel/compilers_and_libraries_2019.0.117/linux/mkl

all: aout

aout: mod_afssh.o AFSSH.o
	ifort -o aout mod_afssh.o AFSSH.o -qopt-matmul -ipo -O3 -no-prec-div -static-intel -fp-model fast=2 -xHost -mkl

%.o: %.f90
	ifort -c $< -qopt-matmul -ipo -O3 -no-prec-div -static-intel -fp-model fast=2 -xHost

quick:
	gfortran -o aout mod_afssh.f90 AFSSH.f90  ~/lapack-3.8.0/liblapack.a ~/lapack-3.8.0/librefblas.a

clean:
	rm *.o aout

