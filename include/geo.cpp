#include <math.h>
#include <complex>
#include "geo.h"
// #include <boost/math/special_functions/hankel.hpp>
#include "mex.h"
#include <sstream>
using namespace std;
const double pi= 4.0*atan(1.0);

typedef complex<double> dcomp;
const dcomp i(0.0,1.0);

const int  Nd = 2;



////////////////Geometria dada la parametrizacion/////////////////////////
void curves_c(double* x, double* y,double t,int curv){
	*x=cos(pi/(double)Nd*(t+2*curv-1));
	*y=sin(pi/(double)Nd*(t+2*curv-1));
}

void curves_r(double* x, double* y,double t,int curv){
	*x=0.5*(t+1)*cos(2*pi/(double)Nd*(curv-1));
	*y=0.5*(t+1)*sin(2*pi/(double)Nd*(curv-1));
}

void tang_curves_c(double* x, double* y,double t,int curv){
	*x=-sin(pi/(double)Nd*(t+2*curv-1))*(pi/(double)Nd);
	*y=cos(pi/(double)Nd*(t+2*curv-1))*(pi/(double)Nd);
}

void tang_curves_r(double* x, double* y,double t,int curv){
	*x=0.5*cos(2*pi/(double)Nd*(curv-1));
	*y=0.5*sin(2*pi/(double)Nd*(curv-1));
}



void tangp(double* x, double* y,double t,int curv){
}


void geomp(double* x, double* y, double t, int domt, int dom_op){

        mxArray *params[3];

    params[0] = mxCreateDoubleScalar(t);

    params[1] = mxCreateDoubleScalar(domt);

    params[2] = mxCreateDoubleScalar(dom_op);

    mxArray *coords[2];

    mexCallMATLAB(2, coords, 3, params, "GeoWrapP");

    *x = mxGetScalar(coords[0]);

    *y = mxGetScalar(coords[1]);
//
      mxDestroyArray(params[0]);

      mxDestroyArray(params[1]);
//
      mxDestroyArray(params[2]);
// //
      mxDestroyArray(coords[0]);
// //
      mxDestroyArray(coords[1]);


}

void geom(double* x, double* y, double t, int domt, int dom_op){
  //entrega la geometria como funcion de t en -1,1.
//     if(dom_op==0){
// 	curves_c( x,  y, -t, domt+1);
//     }else{
//         if(domt==0){
//              curves_r( x,  y, t, dom_op);
//         }else if(domt==1){
//          curves_c( x,  y,t, dom_op);
//         }else if(domt==2){
//              curves_r(x,  y, -t,dom_op+1);
//         }
//     }

    mxArray *params[3];

    params[0] = mxCreateDoubleScalar(t);

    params[1] = mxCreateDoubleScalar(domt);

    params[2] = mxCreateDoubleScalar(dom_op);

    mxArray *coords[2];

    mexCallMATLAB(2, coords, 3, params, "GeoWrap");

    *x = mxGetScalar(coords[0]);

    *y = mxGetScalar(coords[1]);
//
      mxDestroyArray(params[0]);

      mxDestroyArray(params[1]);
//
      mxDestroyArray(params[2]);
// //
      mxDestroyArray(coords[0]);
// //
      mxDestroyArray(coords[1]);

}



double inv_geom(double x, double y, int domt, int dom_op){
}



void geompp(double* x, double* y, double t, int domt, int dom_op){

      mxArray *params[3];

    params[0] = mxCreateDoubleScalar(t);

    params[1] = mxCreateDoubleScalar(domt);

    params[2] = mxCreateDoubleScalar(dom_op);

    mxArray *coords[2];

    mexCallMATLAB(2, coords, 3, params, "GeoWrapPP");

    *x = mxGetScalar(coords[0]);

    *y = mxGetScalar(coords[1]);
//
      mxDestroyArray(params[0]);

      mxDestroyArray(params[1]);
//
      mxDestroyArray(params[2]);
// //
      mxDestroyArray(coords[0]);
// //
      mxDestroyArray(coords[1]);


}

int* trial_test(int dom_op_trial, int dom_op_test,int Dt_trial, int Dt_test){

     int *tt =  new int[Dt_trial];

    double *aux;

    mxArray *params[2];

    params[0] = mxCreateDoubleScalar(dom_op_trial);

    params[1] = mxCreateDoubleScalar(dom_op_test);

    mxArray *lhs[1];

    mexCallMATLAB(1, lhs, 2, params, "TrialTest");

//     for(int ii(0);ii<Dt_trial;++ii)
//     {
//         tt[ii] = (int)mxGetScalar(lhs[ii]);
// //
// //         mxDestroyArray(lhs[ii]);
//     }
//
    aux = mxGetPr(lhs[0]);


    for(int ii(0);ii<Dt_trial;++ii)
    {
        tt[ii] = (int)aux[ii];
//
//         mxDestroyArray(lhs[ii]);
    }

    mxDestroyArray(params[0]);

    mxDestroyArray(params[1]);

    return tt;
}




int reg( int domt, int doms, int dom_op){
    //funcion que determina si se regulariza o no la integral W3, que corresponde
    //a la integral de los valores de frontera (t=-1,t=1) de Vu_trial(t)
    // 2 se regulariza completamente
    // 1/-1 se regulariza solo el -1 o el 1 (entrega cual se regulariza)
    //-2 es disjunta no es necesario regularizar.
    // 0 si se regularian ambos pero no coinsiden los dominios (caso circulo)
     if(domt==doms){
         return 2;
     }else if ((dist(1,-1,domt,doms,dom_op)<pow(10,-10))&&(dist(-1,1,domt,doms,dom_op)<pow(10,-10))){
         return 0;
     }else if(dist(1,-1,domt,doms,dom_op)<pow(10,-10)){
         return 1;
     }else if(dist(-1,1,domt,doms,dom_op)<pow(10,-10)){
         return -1;
     }else{
         return -2;
     }
}

//Dummy function.
// int* trial(int Dt){
//  int *tt =new int[Dt];
//  for(int j=0; j< Dt; j++){
//      tt[j]=j;
//  }
//  return tt;
// }

////////////////////////////////////////////////////////////////////////////


//////Jacobianos///////////////////////////////////////////////////////////

double J(double t,int domt, int dom_op){
 double x,y;
 geomp(&x,&y,t,domt,dom_op);
 return sqrt(pow(x,2)+pow(y,2));
}

double Jp(double t,int domt, int dom_op){

  double x,y;

 geomp(&x,&y,t,domt,dom_op);

 double J= sqrt(pow(x,2)+pow(y,2));

 double xp, yp;

 geompp(&xp,&yp,t,domt,dom_op);

 return (xp*x+yp*y)/J;

}


///////Otras funciones/////////////////////////////////////////////////////

double  dist(double t, double s, int domt, int doms, int dom_op){
    double xt,yt,xs,ys;
    geom(&xt,&yt,t,domt,dom_op);
    geom(&xs,&ys,s,doms,dom_op);
    return sqrt(pow(xt-xs,2)+pow(yt-ys,2));
}



void normal(double* nx, double* ny,double t, int domt, int dom_op){
    double x,y,j;
    geomp(&x,&y,t,domt,dom_op);
    j=sqrt(pow(x,2)+pow(y,2));
//      if(dom_op==0){
//          *nx =-y/j;
//          *ny =x/j;
//      }else{
//          *nx =y/j;
//          *ny =-x/j;
//      }
   *nx =y/j;
   *ny =-x/j;
}
