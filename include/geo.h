#include <complex>
using namespace std;

extern void curves_c(double* x, double* y,double t,int curv);
extern void curves_r(double* x, double* y,double t,int curv);
extern void tang_curves_c(double* x, double* y,double t,int curv);
extern void tang_curves_r(double* x, double* y,double t,int curv);
extern void geom(double* x, double* y, double t, int domt, int dom_op);
extern void geomp(double* x, double* y, double t, int domt, int dom_op);
extern void geompp(double* x, double* y, double t, int domt, int dom_op);
extern int reg(int domt, int doms, int dom_op);
extern int* trial_test(int dom_op_trial, int dom_op_test,int Dt_trial, int Dt_test); 
// extern int* trial(int Dt);
extern double J(double t,int domt, int dom_op);
extern double Jp(double t,int domt, int dom_op);
extern double  dist(double t, double s, int domt, int doms, int dom_op);
extern void normal(double* nx, double* ny,double t, int domt, int dom_op);
extern double inv_geom(double x, double y, int domt, int dom_op);
extern void tangp(double* x, double* y,double t,int curv);