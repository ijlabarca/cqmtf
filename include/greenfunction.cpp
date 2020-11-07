// #include <boost/math/special_functions/hankel.hpp>
#include <complex_bessel.h>
#include "greenfunction.h"
#include <vector>


const double pii= 4.0*atan(1.0);

dcomp FreeSpace::IntegralKernel()
{
    // return 0.25*i*boost::math::cyl_hankel_1(0,m_WaveNumber.real()*m_distance);
    return 0.25*i*sp_bessel::hankelH1(0, i*m_WaveNumber*m_distance );
}

dcomp GreenFunctionBase::evaluateRegularizator(double t, double s)
{

    return dcomp(-1/(2*pii)*log(abs(t-s)),0);

}


dcomp FreeSpace::evaluateGF( int n,int typen, int l,int typel,
        int domt, int doms,int dom_op)
{
    SetPointsFromGeo(n,typen,l,typel,domt,doms,dom_op);

    ComputeDistance();

    return IntegralKernel();

}


 dcomp FreeSpace::IntegralKernelDev()
 {
     // return -m_WaveNumber.real()*0.25*i*boost::math::cyl_hankel_1
             // (1,m_WaveNumber.real()*m_distance);

     return m_WaveNumber*0.25*sp_bessel::hankelH1
              (1, i*m_WaveNumber*m_distance );

 }


 void FreeSpace::ComputeGradDistance()    // (x-y)/||x-y||
 {
    m_gradDistance[0] = (m_x[0]-m_y[0])/m_distance;

    m_gradDistance[1] = (m_x[1]-m_y[1])/m_distance;

 }

 dcomp FreeSpace::evaluateGFDetNx(int n,int typen, int l,int typel,
         int domt, int doms,int dom_op)
 {
    SetPointsFromGeo(n,typen,l,typel,domt,doms,dom_op);

    ComputeDistance();

    ComputeGradDistance();

    SetXNormalFromGeo(n,typen,domt,dom_op);

    double NormaldotGrad = m_xNormal[0]*m_gradDistance[0]+
            m_xNormal[1]*m_gradDistance[1];

    return IntegralKernelDev()*NormaldotGrad;

 }

  dcomp FreeSpace::evaluateGFDetNy(int n,int typen, int l,int typel,
         int domt, int doms,int dom_op)
 {
    SetPointsFromGeo(n,typen,l,typel,domt,doms,dom_op);

    ComputeDistance();

    ComputeGradDistance();

    SetYNormalFromGeo(l,typel,doms,dom_op);

    double NormaldotGrad = m_yNormal[0]*(-m_gradDistance[0])+
            m_yNormal[1]*(-m_gradDistance[1]);

    return IntegralKernelDev()*NormaldotGrad;

 }
