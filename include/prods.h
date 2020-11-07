#include <complex>
#include "domain.h"

using namespace std;
typedef complex<double> dcomp;
extern dcomp pV(int m, int l ,dcomp** ft_cofs,int Ng,double* xg,double* wg, dcomp k,int domt,int doms);
extern dcomp pK(int m, int l ,dcomp** ft_cofs,int Ng,double* xg,double* wg, dcomp k,int domt,int doms);
extern dcomp pW(int m, int l ,
        Domain & dom,
        dcomp*** ft_cofs,dcomp** wb,int Ng,double* xg,double* wg, dcomp k,int domt,int doms,int dom_op);
//extern double p1(int m, int l,int Nc,double xc[],double wc[], int domt, int dom_op);
extern double p1(dcomp* ft_cofs, int l);
extern dcomp pld(int l,dcomp* ft_cofs,int dom_op);
extern dcomp pld_d(int l,int Nc,double xc[],double wc[], int doms, int dom_op,dcomp k,int ll, double R);
extern dcomp pld_n(int l,int Nc,double xc[],double wc[], int doms, int dom_op,dcomp k, int ll, double R);
extern double *t_polynomial ( int m, int n, double x[] );
extern double chebT ( int n, double x );
extern double *u_polynomial ( int m, int n, double x[] );
extern double chebU ( int n, double x );
extern double chebU_p ( int n, double x );
