#include <complex>
#include "mex.h"
#include "domain.h"
#include "IncidentField.h"

using namespace std;
typedef complex<double> dcomp;

extern void Adom(dcomp** A, dcomp k, int* N,int Nc, int Ng, double*xg, double* wg,int dom_op, int Dt, int NT, Domain &dom);

extern void Xdoms(double** X,  int* N_trial,int* N_test,int dom_op_trial,
        int Dt_trial,int dom_op_test, int Dt_test,
        int NT_trial, int NT_test,
        Domain &domTrial,mxArray *doms[2]);

extern void bdom(dcomp* b, int* N,int dom_op, int Dt,int Dt_0,
        IncidentField& inc,int NT,Domain &dom0,mxArray *doms[2]);

extern void  bdomPlaneWave(dcomp* b, int* N,int dom_op, int Dt,int Dt_0,
        dcomp k,double alpha,int NT,Domain &dom0,mxArray *doms[2]);

extern int* trial_test(int Dt_trial,mxArray *doms[2]);

int* trial(int Dt);
