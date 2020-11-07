#include <complex>
#include <fftw3.h>
using namespace std;
typedef complex<double> dcomp;

/*
-----------------------------------------------------------
-----------------------------------------------------------
Tchv_cof
Computes the Chebyschev coefficients associated to fxc.
-----------------------------------------------------------
n       : length of fxc
fxc     : std::complex<double> pointer with values
          of f evaluated at cheb.points xc
p       : FFTW Plan
-----------------------------------------------------------
-----------------------------------------------------------
*/
extern void *Tchv_cof(int n, dcomp* fxc,fftw_plan p);
extern fftw_plan plan_0(int n, dcomp* fxc);
extern  double *chebpts(int n);
