#ifndef MY_HEADER_H2
#define MY_HEADER_H2

#include "mex.h"
#include <math.h>

using namespace std;

class Domain {

protected:
    int m_NumIntrfz;
    int m_NGauss;
    int m_NCheb;
    double *m_Gauss;
    double *m_Cheb;

public:
        Domain( int NIntrz,
            int Ng,
            int Nc,
            double *Gauss,
            double *Cheb):
                m_NumIntrfz(NIntrz),
                m_NGauss(Ng),
                m_NCheb(Nc),
                m_Gauss(Gauss),
                m_Cheb(Cheb)

    {
    }

     ~Domain(){}

    void Geo(double *x, double *y,int n,int type,int intfs);

    double Jac(int n,int type,int intfs);

    double JacP(int n,int type,int intfs);

    void Normal(double *x, double *y,int n,int type,int intfs);

    double dist(int n,int typen, int l,int typel, int intfst, int intfss);

    int reg(int intfst, int intfss);

    int GetNc(){return m_NCheb;}

    int GetNg(){return m_NGauss;}

    int GetNItrfz(){return m_NumIntrfz;}

};
 #endif
