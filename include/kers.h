#include <complex>
#include "greenfunction.h"
#include "domain.h"

using namespace std;
typedef complex<double> dcomp;
extern void pk(double* dk, double* pk, double t, double s,  int domt, int doms, int dom_op);
extern void pka(double* dk, double* pk, double t, double s,  int domt, int doms, int dom_op);
extern void pw(double* dw, double* pw, double t, double s,  int domt, int doms, int dom_op);





extern dcomp G0(double t, double s);

extern dcomp GB(double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf);

extern dcomp KB(int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf);

extern dcomp KAB(int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf);

extern dcomp WB(double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf,
        Domain & dom);

extern dcomp GBB(double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf);

extern dcomp WA1( double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf,
        Domain & dom);

extern dcomp WA2( double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf,
        Domain & dom);

extern dcomp WA3( double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf,
        Domain & dom);

extern dcomp WA4( double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf,
        Domain & dom);

extern dcomp WB1( int sig,double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf,
        Domain & dom);

extern dcomp WB2(int sig, double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf,
        Domain & dom);


extern dcomp uinc_d(dcomp k0, int l,int type, int doms ,int dom_op,
        int ll, double R,
        Domain & dom0);

extern dcomp uinc_n(dcomp k0, int l,int type, int doms ,int dom_op,
        int ll, double R,
        Domain & dom0);

extern dcomp planeWave(dcomp k0, int l,int type, int doms ,double alpha,
        Domain & dom0);

extern dcomp planeWaveN(dcomp k0, int l,int type, int doms ,double alpha,
        Domain & dom0);

extern double Argument(double x, double y);
