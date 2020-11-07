#include "domain.h"

//eliminar los ifs y ahcer 2 funciones para que sea mas rapido.


/*
n       : interface
type    : 0: Gauss,   1: Chebyschev

m_Gauss = [x, y, nx, ny, Jac, JacP] x interface

m_Cheb = [x, y, nx, ny, Jac, JacP] x interface
*/


void Domain::Geo(double *x, double *y,int n,int type,int intfs)
{


    if (type ==0)
        //gauss
    {
        *x = m_Gauss[intfs+6*m_NumIntrfz*n];

        *y = m_Gauss[intfs+m_NumIntrfz+6*m_NumIntrfz*n];

    }

    else

    {
        *x = m_Cheb[intfs+6*m_NumIntrfz*n];

        *y = m_Cheb[intfs+m_NumIntrfz+6*m_NumIntrfz*n];
    }


}

double Domain::Jac(int n,int type,int intfs)
{

    if (type ==0)
    {
        return m_Gauss[intfs+4*m_NumIntrfz+6*m_NumIntrfz*n];

    }
    else
    {
        return m_Cheb[intfs+4*m_NumIntrfz+6*m_NumIntrfz*n];

    }

}

double Domain::JacP(int n,int type,int intfs)
{

    if (type ==0)
    {
        return m_Gauss[intfs+5*m_NumIntrfz+6*m_NumIntrfz*n];

    }
    else
    {
        return m_Cheb[intfs+5*m_NumIntrfz+6*m_NumIntrfz*n];

    }

}

void Domain::Normal(double *x, double *y,int n,int type,int intfs)
{

    if (type ==0)
        //gauss
    {
        *x = m_Gauss[intfs+2*m_NumIntrfz+6*m_NumIntrfz*n];

        *y = m_Gauss[intfs+3*m_NumIntrfz+6*m_NumIntrfz*n];

    }

    else

    {
        *x = m_Cheb[intfs+2*m_NumIntrfz+6*m_NumIntrfz*n];

        *y = m_Cheb[intfs+3*m_NumIntrfz+6*m_NumIntrfz*n];
    }


}

double Domain::dist(int n,int typen, int l,int typel, int intfst, int intfss)
{
    double xt,yt,xs,ys;

    Geo(&xt,&yt,n,typen,intfst);

    Geo(&xs,&ys,l,typel,intfss);

    return sqrt(pow(xt-xs,2)+pow(yt-ys,2));

}

int Domain::reg(int intfst, int intfss)
{
    //funcion que determina si se regulariza o no la integral W3, que corresponde
    //a la integral de los valores de frontera (t=-1,t=1) de Vu_trial(t)
    // 2 se regulariza completamente
    // 1/-1 se regulariza solo el -1 o el 1 (entrega cual se regulariza)
    //-2 es disjunta no es necesario regularizar.
    // 0 si se regularian ambos pero no coinsiden los dominios (caso circulo)

     if(intfst==intfss){
         return 2;
     }
     else if ((dist(m_NCheb-1,1,0,1,intfst,intfss)<pow(10,-10))
     &&(dist(0,1,m_NCheb-1,1,intfst,intfss)<pow(10,-10))){
         return 0;
     }
     else if(dist(m_NCheb-1,1,0,1,intfst,intfss)<pow(10,-10)){
         return 1;
     }
     else if(dist(0,1,m_NCheb-1,1,intfst,intfss)<pow(10,-10)){
         return -1;
     }
     else{
         return -2;
     }

}
