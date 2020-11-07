#include <math.h>
#include <complex>
#include "kers.h"
#include "geo.h"
#include "domain.h"
// #include <boost/math/special_functions/hankel.hpp>
#include <complex_bessel.h>
#include "mex.h"


#include "greenfunction.h"

using namespace std;
const double pi= 4.0*atan(1.0);

///////////funciones auxiliares (multiplicadores de kernells)/////////////////

void pk(double* dk, double* pk, double t, double s,  int domt, int doms, int dom_op){
    double xt,yt,xs,ys,nx,ny;
    geom(&xt,&yt,t,domt,dom_op);
    geom(&xs,&ys,s,doms,dom_op);
    normal(&nx,&ny,t,domt,dom_op);
    *dk=sqrt(pow(xt-xs,2)+pow(yt-ys,2));
    *pk=((xt-xs)*nx+(yt-ys)*ny)/(*dk);
}

void pka(double* dk, double* pk, double t, double s,  int domt, int doms, int dom_op){
    double xt,yt,xs,ys,nx,ny;
    geom(&xt,&yt,t,domt,dom_op);
    geom(&xs,&ys,s,doms,dom_op);
    normal(&nx,&ny,s,doms,dom_op);
    *dk=sqrt(pow(xt-xs,2)+pow(yt-ys,2));
    *pk=((xs-xt)*nx+(ys-yt)*ny)/(*dk);
}

void pw(double* dw, double* pw, double t, double s,  int domt, int doms, int dom_op){
    double nxt,nyt,nxs,nys;
    normal(&nxt,&nyt,t,domt,dom_op);
    normal(&nxs,&nys,s,doms,dom_op);
    *dw=dist(t,s,domt,doms,dom_op);
    *pw=nxt*nxs+nyt*nys;
}


////////////////////////////////////////////////////////////////////////////////////////

//////Kernels//////////////////////////////////////////////////////////////////


dcomp G0(double t, double s){
	return dcomp(-1/(2*pi)*log(abs(t-s)),0);
}

dcomp GB(double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf){
	if(domt==doms){
		return gf.evaluateGF(n,typen,l,typel,domt,doms,dom_op)
        - gf.evaluateRegularizator(t,s);
	}else{
		return gf.evaluateGF(n,typen,l,typel,domt,doms,dom_op);
	}
}

dcomp KB(int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf){

    return gf.evaluateGFDetNx(n,typen,l,typel,domt,doms,dom_op);
}

dcomp KAB(int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf){

    return gf.evaluateGFDetNy(n,typen,l,typel,domt,doms,dom_op);
}

dcomp WB(double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf,
        Domain & dom){

   	double p;

    double nx1,ny1,nx2,ny2;

    dom.Normal(&nx1,&ny1,n,typen,domt);

    dom.Normal(&nx2,&ny2,l,typel,doms);

    p = nx1*nx2+ny1*ny2;

    dcomp k = gf.GetWaveNumber();

	if(domt==doms){
		return -k*k*
                (-1.0*gf.evaluateGF(n,typen,l,typel,domt,doms,dom_op)*p
                +gf.evaluateRegularizator(t,s));
	}else{
		return k*k*
                gf.evaluateGF(n,typen,l,typel,domt,doms,dom_op)*p;
	}
}

dcomp GBB(double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf){

	return gf.evaluateGF(n,typen,l,typel,domt,doms,dom_op)
        - gf.evaluateRegularizator(t,s);

}

dcomp WA1( double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf,
        Domain & dom){

    double jt,js,p;

    jt= dom.Jac(n,typen,domt);

    js= dom.Jac(l,typel,doms);

    p=1.0/(jt*js);

	if(domt==doms){

        double p0 = pow(jt,-2);

		return gf.evaluateGF(n,typen,l,typel,domt,doms,dom_op)*p
                -gf.evaluateRegularizator(t,s)*p0;
	}else{
		return gf.evaluateGF(n,typen,l,typel,domt,doms,dom_op)*p;
	}
}

dcomp WA2( double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf,
        Domain & dom){

    double js;

    js= dom.Jac(l,typel,doms);

	if(domt==doms){

        double jt =dom.Jac(n,typen,domt);

		return gf.evaluateGF(n,typen,l,typel,domt,doms,dom_op)/js-
               gf.evaluateRegularizator(t,s)/jt;
	}else{
		return gf.evaluateGF(n,typen,l,typel,domt,doms,dom_op)/js;
	}
}

dcomp WA3( double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf,
        Domain & dom){
    double js,jsp,p;

    js= dom.Jac(l,typel,doms);

    jsp=dom.JacP(l,typel,doms);

    p=jsp*pow(js,-2);

	if(domt==doms){

        double jt,jtp,p0;

        jt=dom.Jac(n,typen,domt);

        jtp=dom.JacP(n,typen,domt);

        p0 = jtp*pow(jt,-2);

		return gf.evaluateGF(n,typen,l,typel,domt,doms,dom_op)*p
                -gf.evaluateRegularizator(t,s)*p0;

	}else{

		return gf.evaluateGF(n,typen,l,typel,domt,doms,dom_op)*p;
	}
}

dcomp WB1(int sig, double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf,
        Domain & dom){


    double js;
    int dmx;
    int sign=1;
    int bol=dom.reg(domt,doms);

    int ns = n;

    if((abs(bol)==1)&&(sig*bol==-1)){
        bol=-2;
    }
    if(bol>-2){
        dmx=doms;
    }else{
        dmx=domt;
    }
    if(abs(bol)<2){
        sign=-1;
    }

    if (sign == -1)
    {
        ns = dom.GetNg()-1-n;
    }

    js=dom.Jac(l,typel,doms);

//      if( (dom_op==0)&&(domt==0)&&(doms==1))
//     {
//           double jt=dom.Jac(ns,typen,dmx);
//
//         dcomp aux =gf.evaluateGF(ns,typen,l,typel,dmx,doms,dom_op)/js;
//
//         mexPrintf("Sign %d domx %d doms %d t %f s %f | Gf: %f , %f  \n",sign, dmx,doms,t,s,aux.real(),aux.imag());
//     }

    if(bol>-2){

        double jt=dom.Jac(ns,typen,dmx);

        return gf.evaluateGF(ns,typen,l,typel,dmx,doms,dom_op)/js-
                gf.evaluateRegularizator((double)sign*t,s)/jt;

    }else{
        return gf.evaluateGF(ns,typen,l,typel,dmx,doms,dom_op)/js;
    }

}



dcomp WB2(int sig, double t, double s,
        int n,int l,int typen,int typel,
        int domt, int doms, int dom_op,
        GreenFunctionBase & gf,
        Domain & dom){

    double js,jsp,p;
    int dmx,ns;

    ns =n;

    int sign=1;
    int bol=dom.reg(domt,doms);
    if((abs(bol)==1)&&(sig*bol==-1)){
        bol=-2;
    }
    if(bol>-2){
        dmx=doms;
    }else{
        dmx=domt;
    }
    if(abs(bol)<2){
        sign=-1;
    }

    if (sign == -1)
    {
        ns = dom.GetNg()-1-n;
    }

    js=dom.Jac(l,typel,doms);

    jsp=dom.JacP(l,typel,doms);

    p=jsp*pow(js,-2);

    if(bol>-2){

        double jt,jtp,p0;

        jt=dom.Jac(ns,typen,dmx);

        jtp=dom.JacP(ns,typen,dmx);

        p0=jtp*pow(jt,-2);


        return gf.evaluateGF(ns,typen,l,typel,dmx,doms,dom_op)*p-
                      gf.evaluateRegularizator((double)sign*t,s)*p0;

    }else{

        return gf.evaluateGF(ns,typen,l,typel,dmx,doms,dom_op)*p;

    }
}

// k complejo
dcomp uinc_d(dcomp k0, int l,int type, int doms ,int dom_op,
        int ll, double R,
        Domain &dom0){


    double x,y,norm,unitx,unity,argx;

    dom0.Geo(&x,&y,l,type,doms);

//     geom(&x,&y,s,doms,0);
//
    norm = sqrt(x*x+y*y);

    unitx=x/norm;

    unity=y/norm;

     argx = Argument(unitx,unity);

//     argx = atan(unity/unitx);

    dcomp H,Y;

    dcomp J;

//     H=boost::math::cyl_hankel_1(abs(l),k0*R);

    H=1.0;

    // J= boost::math::cyl_bessel_j(ll,k0*norm);
    J = sp_bessel::besselJ(ll,i*k0*norm);

    Y = 1/sqrt(2*pi)*exp(i*(double)ll*argx);

    return Y*J;

}

// k complejo
dcomp uinc_n(dcomp k0, int l,int type, int doms ,int dom_op,
        int ll, double R,
        Domain &dom0){

    double nx,ny;

    dom0.Normal(&nx,&ny,l,type,doms);

    double x,y,norm,unitx,unity,argx;

    dom0.Geo(&x,&y,l,type,doms);

    norm = sqrt(x*x+y*y);

    unitx=x/norm;

    unity=y/norm;

    argx = Argument(unitx,unity);

    dcomp H,Y;

    dcomp J,Jp;

//     H=boost::math::cyl_hankel_1(abs(l),k0*R);

    H=1;

    // J= boost::math::cyl_bessel_j(ll,k0*norm);
    J = sp_bessel::besselJ(ll,i*k0*norm);


    Y = 1.0/sqrt(2*pi)*exp(i*(double)ll*argx);

    if( abs(l)> 0)
        // Jp= 0.5*(boost::math::cyl_bessel_j(ll-1,k0*norm)-
                 // boost::math::cyl_bessel_j(ll+1,k0*norm) );

        Jp= 0.5*(sp_bessel::besselJ(ll-1,i*k0*norm)-
                sp_bessel::besselJ(ll+1,i*k0*norm) );




    else
        // Jp= -1.0*boost::math::cyl_bessel_j(ll+1,k0*norm) ;
        Jp = -1.0*sp_bessel::besselJ(ll+1,i*k0*norm);


    // return (k0*Jp*Y*(unitx*nx+unity*ny)+
            // i*(double)ll/norm*J*Y*(-unity*nx+unitx*ny));


    return (i*k0*Jp*Y*(unitx*nx+unity*ny)+
            i*(double)ll/norm*J*Y*(-unity*nx+unitx*ny));
   }
///////////////////////////////////////////////////////////////////////

// k complejo
dcomp planeWave(dcomp k0, int l,int type, int doms ,double alpha,
        Domain &dom0){

    double x,y;
    dom0.Geo(&x,&y,l,type,doms);
    // return exp(i*k0*(x*cos(alpha)+y*sin(alpha)));
    return exp(-1.0*k0*(x*cos(alpha)+y*sin(alpha)));


}

// k complejo
dcomp planeWaveN(dcomp k0, int l,int type, int doms ,double alpha,
        Domain &dom0){


    double nx,ny;
    dom0.Normal(&nx,&ny,l,type,doms);
    // return i*k0*planeWave(k0,l,type,doms,alpha,dom0)*(cos(alpha)*nx+sin(alpha)*ny);
    return -1.0*k0*planeWave(k0,l,type,doms,alpha,dom0)*(cos(alpha)*nx+sin(alpha)*ny);


}


double Argument(double x, double y)
{

    if((x>=0)&&(y>=0))
        return atan(y/x);
    else if((x<0)&&(y>0))
        return pi-atan(y/(-x));
    else if((x>0) &&(y<0))
        return 2*pi-atan((-y)/x);
    else
        return pi+atan(y/x);



}
