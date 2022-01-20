FC=mpif90 
FFLAGS=-O3 -Wall -Wextra -fdefault-real-8 -fimplicit-none -std=f2003 
INCLUDE=-I/usr/include -I./2decomp_fft/include
LDFLAGS=-L./2decomp_fft/lib -lm -lfftw3 -l2decomp_fft 

all: main

main.o: main.f90
	$(FC) $(INCLUDE) $(FFLAGS) -c main.f90

main: main.o
	$(FC) -o main main.o $(LDFLAGS) -lfftw3

clean:
	rm -f *.o *.mod main *.tec *.dat
