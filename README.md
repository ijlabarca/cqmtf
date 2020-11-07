# CQ - MTF
Convolution Quadrature and Multiple traces formulation using Galerkin discretization with spectral elements

MATLAB files are provided. Computations are done with C++ and linked to MATLAB through a MEX file. 

|C++ files:           | Description |
| ------------------- | ---------------|
| `mex_funs.cpp`      |  MATLAB link, calls A or X |
| `Adom.cpp`          |    |
| `prods.cpp`         |   |
| `kers.cpp`          |  |
| `greenfunction.cpp` | Kernel of BIOs |


|MATLAB files:           | Description |
| ------------------- | ---------------|
| `solverMtf28.m`      |  Class for solver |
| `main_CQMTF.m`       | Computes solution for a problem with BDF2-CQ, plots in time domain   |
| `main_RKCQ_MTF.m`    | Computes solution for a problem with RKCQ, plots in time domain   |
|         |   |
|           |  |
|  |  |

The following libraries are required for the proper use of the code:

### complex_bessel
A C++ library to evaluate Bessel functions of all kinds. More information can be found on the website.
https://github.com/joeydumont/complex_bessel

### FFTW
http://www.fftw.org/


