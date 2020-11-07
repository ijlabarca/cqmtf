using namespace std;
#include "domain.h"
#include "IncidentField.h"

typedef complex<double> dcomp;
extern void in_cells(dcomp k,dcomp***** Opd,dcomp**** Wb,int Dt, int dom_op, int Ng, double*xg, int Nc,double* xc,double xp, Domain &dom);
extern void in_cells0(dcomp*** uin,int Dt,dcomp k, double alpha, int Nc,double* xc);
extern void in_cells_d(dcomp*** Opd,int Dt, int dom_op_trial,int dom_op_test,int Nc,double* xc,int* N,Domain &domTrial);
extern void in_cells_ld(dcomp*** Opd,int Dt, int Nc,double* xc,IncidentField& inc,Domain&dom0);

extern void  in_cells_ldPlaneWave(dcomp*** Opd,int Dt, int Nc,double* xc,dcomp k, double alpha,Domain &dom0);
