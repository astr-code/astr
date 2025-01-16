# ASTR code 
Version 2.0 

ASTR code is a high-order finite-difference flow solver for compressible turbulence research. This project explores the usage of CUDA-Fortran to parallelise the ASTR code. The tgv solver for a 3D case has been parallelised using CUDA.

# Download, Installation and Compilation
Required dependencies: Fortran 90, NVIDIA HPC SDK, CMAKE

The installation guide for NVIDIA HPC SDK can be found at [Installation Guide](https://docs.nvidia.com/hpc-sdk/hpc-sdk-install-guide/index.html)

## Download the gpu accelerated astr code:

Clone the git repository:
```
$ git clone https://github.com/terencel411/astr.git
```
## Compilation:
The `Makefile` gives a complete and safe way of compiling the code.

Go to the directory where the miniapps code is present

```
$ cd astr/miniapps/tgv_solver_3d
```

The cpu and the gpu accelerated codes (`tgvsolver_cpu.f90` and `tgvsolver_gpu.cuf`) are present in the same directory, which can be compiled using the following cmake commands

Compile and execute the cpu & gpu code

```
$ cmake all
```

The cpu and gpu code can also be compiled and executed separately

```
$ cmake cpu
$ cmake gpu
```

## Acceleration Comparison
The time acceleration statistics can be obtained by running the following command. 

```
$ cmake compare
```

A text file `time_report.txt` will be generated with the accelerations statistics.






