#ifndef MY_HEADER_H
#define MY_HEADER_H
#include "mex.h"

// #include "geo.h"
#include "domain.h"
#include <vector>

using namespace std;

typedef complex<double> dcomp;

const dcomp i(0.0,1.0);

class GreenFunctionBase {

protected:
    dcomp m_WaveNumber;
    vector<double> m_x;
    vector<double> m_y;
    vector<double> m_xNormal;
    vector<double> m_yNormal;
    double m_distance;
    Domain* m_dom;

public:

    // Recordar que k corresponde a Δu - k^2 u = 0, k complejo
    GreenFunctionBase(dcomp k, Domain *dom):
    m_WaveNumber(k),
    m_dom(dom)
    {
        m_x.resize(2,0.0);
        m_y.resize(2,0.0);
        m_xNormal.resize(2,0.0);
        m_yNormal.resize(2,0.0);
    }

    void setTwoDVector(double v1, double v2, vector<double>& vec){
      vec[0] = v1;  vec[1] = v2;}

    void setXpoint  (double x1,double x2){setTwoDVector(x1,x2,m_x);}

    void setYpoint  (double y1,double y2){setTwoDVector(y1,y2,m_y);}

    void setXNormal  (double x1,double x2){setTwoDVector(x1,x2,m_xNormal);}

    void setYNormal  (double x1,double x2){setTwoDVector(x1,x2,m_yNormal);}

    void setXpoint  (vector<double> data ){m_x = data;}

    void setYpoint  (vector<double> data){m_y = data;}

    void setXNormal  (vector<double> data){m_xNormal = data;}

    void setYNormal  (vector<double> data){m_yNormal = data;}

    void SetXFromGeo(int n,int type, int domt, int dom_op){
        double x1,x2;
        m_dom->Geo(&x1, &x2, n,type,domt);
        setXpoint(x1,x2);
    }

    void SetYFromGeo(int  n,int type, int doms, int dom_op){
        double y1,y2;
        m_dom->Geo(&y1, &y2, n,type,doms);
        setYpoint(y1,y2);
    }

    void SetXNormalFromGeo(int n,int type, int domt, int dom_op)
    {
        double nx1,nx2;
        m_dom->Normal(&nx1, &nx2, n,type,domt);
        setXNormal(nx1,nx2);
    }

    void SetYNormalFromGeo(int n,int type, int doms, int dom_op)
    {
        double ny1,ny2;
        m_dom->Normal(&ny1, &ny2, n,type,doms);
        setYNormal(ny1,ny2);
    }

    void ComputeDistance()
    {
        m_distance = sqrt(pow(m_x[0]-m_y[0],2)+
                pow(m_x[1]-m_y[1],2));

        if (m_distance < 0.00000001){
          m_distance = 0.00000001;
        }

    }

    void SetPoints(double x1,double x2,double y1,double y2){
        setXpoint(x1,x2);
        setYpoint(y1,y2);
    }

    void SetPoints(vector<double> datax,vector<double> datay){
        setXpoint(datax);
        setYpoint(datay);
    }

    void SetPointsFromGeo(int n,int typen, int l,int typel,
            int domt, int doms,int dom_op)
    {
        SetXFromGeo(n,typen, domt, dom_op);

        SetYFromGeo(l,typel, doms, dom_op);

    }

    dcomp GetWaveNumber(){return m_WaveNumber;}

     virtual dcomp evaluateRegularizator(double t, double s);

     virtual dcomp evaluateGF( int n,int typen, int l,int typel,
             int domt, int doms,int dom_op)=0;

     virtual dcomp evaluateGFDetNx(int n,int typen, int l,int typel,
             int domt, int doms,int dom_op)=0;

     virtual dcomp evaluateGFDetNy(int n,int typen, int l,int typel,
             int domt, int doms,int dom_op)=0;

};

class FreeSpace : public GreenFunctionBase
{
    protected:

        vector<double> m_gradDistance;

        dcomp IntegralKernel();

        dcomp IntegralKernelDev();

        void ComputeGradDistance();

    public:

        FreeSpace(dcomp k, Domain *dom)
            :GreenFunctionBase(k,dom)
            {
                m_gradDistance.resize(2,0.0);
            }

     virtual dcomp evaluateGF( int n,int typen, int l,int typel,
             int domt, int doms,int dom_op);

     virtual dcomp evaluateGFDetNx(int n,int typen, int l,int typel,
             int domt, int doms,int dom_op);

     virtual dcomp evaluateGFDetNy(int n,int typen, int l,int typel,
             int domt, int doms,int dom_op);

};

#endif
