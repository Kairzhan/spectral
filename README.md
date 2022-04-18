# Spectral code for 3D periodic flow simulations

To succesfully compile the code one needs to install **MPI** (for example https://www.open-mpi.org/), **FFTW3** (http://fftw.org/) and **2decomp&FFT** (http://www.2decomp.org/) libraries. Optional CUDA-CUFFT version requires **NVIDIA HPC SDK** to be installed (https://developer.nvidia.com/hpc-sdk).

## Building

This repository contains two versions of the code. The first, CPU version of the code is stored in `src/` directory. The second, GPU (CUDA) version of the code is stored in `src/cuda/` directory.

To build the CPU version of the code one needs to install/compile all the necessary dependencies first. Under Ubuntu use the following commands to install OpenMPI and FFTW3

```
sudo apt install libfftw3-3 libopenmpi3 openmpi-bin libopenmpi-dev
```

Installation of **2decomp&FFT** is straightforward:

```
wget http://www.2decomp.org/download/2decomp_fft-1.5.847.tar.gz
tar -xzf 2decomp_fft-1.5.847.tar.gz
cd 2decomp_fft
cp src/Makefile.inc.x86 src/Makefile.inc
make
```

To finish the build of spectral solver one needs to fix Makefile with correct path to 2decomp&FFT library. Final step to compile the code is also simple:
```
make clean
make
```

Produced executable file **main** can be copied to work directory. Simulation is started through MPI wrapper:
```
mpirun -np 4 ./main
```

The CUDA version of the code is available under `src/cuda`. To compile the CUDA version of the code first one needs to fix PATH to CUDA libraries. To perform this modify **LDFLAGS** variable in `src/cuda/Makefile`. Then simply do:
```
make clean
make
```
Simulation is started explicitly using binary executable:
```
./main
```

## Usage

`utils` directory contains Python script for postprocessing. During the simulation the code produces 3D binary dumps in `uuu*.dat` `vvv*.dat` `www*.dat` files. These files can be processed by `process.py` script. Example usage of this script:
```
python3 ./process.py -n 64 -s 100 -e 2000 -i 25 -m 1 1 1
```
This script should be invoked from the directory which contains `uuu*` `vvv*` `www*` files. In previous example it is explicitly defined that
mesh resolution is `64`, output is started at timestep `100` (-s 100) and ends at `2000` (-e 2000), output frequency is `25` (-i 25). Such **start**, **end** and **frequency** parameters assume that output contains following files:
```
uuu0000100.dat
uuu0000125.dat
uuu0000150.dat
...
uuu0002000.dat
```

More parameters can be found by investigating the script source.
