# Spectral code for 3D periodic flow simulations

To succesfully compile the code one needs to install **MPI** (for example https://www.open-mpi.org/), **FFTW3** (http://fftw.org/) and **2decomp&FFT** (http://www.2decomp.org/) libraries. Optional CUDA-CUFFT version requires **NVIDIA HPC SDK** to be installed (https://developer.nvidia.com/hpc-sdk).

After fixing Makefile with correct path to FFTW/MPI/2decomp&FFT use the following commands to compile the code:
```
make clean
make
```

The CUDA version of the code is available under `src/cuda`.
