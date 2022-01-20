FC=nvfortran
FFLAGS=-O3 
INCLUDE=
LDFLAGS=-L /opt/nvidia/hpc_sdk/Linux_x86_64/2022/math_libs/lib64/ -lcufft 

all: main

main: main.cuf cufft.cuf
	$(FC) $(LDFLAGS) cufft.cuf main.cuf -lcufft

clean:
	rm -f *.o *.mod main *.tec *.dat
