# Spectral code for 3D periodic flow simulations

To succesfully compile the code one needs to install **MPI** (for example https://www.open-mpi.org/), **FFTW3** (http://fftw.org/) and **2decomp&FFT** (http://www.2decomp.org/) libraries. Optional CUDA-CUFFT version requires **NVIDIA HPC SDK** to be installed (https://developer.nvidia.com/hpc-sdk).

After fixing Makefile with correct path to FFTW/MPI/2decomp&FFT use the following commands to compile the code:
```
make clean
make
```

The CUDA version of the code is available under `src/cuda`.

`utils` directory contains Python script for postprocessing. Example usage of this script:
```
python3 ./process.py -n 64 -s 100 -e 2000 -i 25 -m 1 1 1
```
This script should be invoked from the directory which contains `uuu*` `vvv*` `www*` files. In previous example it is explicitly defined that
mesh resolution is `64`, output is started at timestep `100` and ends at `2000`, output frequency is `25`. More parameters can be found by investigating the script source.
