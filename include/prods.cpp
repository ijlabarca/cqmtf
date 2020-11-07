#include <math.h>
#include <complex>
// #include <boost/math/special_functions/bessel.hpp>
#include "kers.h"
#include "prods.h"
#include "geo.h"
#include "domain.h"
using namespace std;
const double pi= 4.0*atan(1.0);


typedef complex<double> dcomp;
//const dcomp i(0.0,1.0);


// k complejo
dcomp pV(int m, int l ,dcomp** ft_cofs,int Ng,double* xg,double* wg, dcomp k,int domt,int doms){
 dcomp u(0,0);
 if(l==0){
    for(int j=0; j<Ng;j++){u+=0.5*pi*wg[j]*(ft_cofs[j][0]-0.5*ft_cofs[j][2])*chebU(m,xg[j]);}
 }else{
    for(int j=0; j<Ng;j++){u+=0.25*pi*wg[j]*(ft_cofs[j][l]-ft_cofs[j][l+2])*chebU(m,xg[j]);}
 }
 if(domt==doms){
     if(l==0){
        for(int j=0; j<Ng;j++){u+=0.25*wg[j]*(log(2)-chebT(2,xg[j])/2)*chebU(m,xg[j]);}
     }else{
        for(int j=0; j<Ng;j++){u+=0.25*wg[j]*(chebT(l,xg[j])/((double)l)-chebT(l+2,xg[j])/((double)(l+2)))*chebU(m,xg[j]);}
     }
 }
 return u;
}

// k complejo
dcomp pK(int m, int l ,dcomp** ft_cofs,int Ng, double* xg,double* wg, dcomp k,int domt,int doms){
 //for k and k' since the trial basis are the same.
 dcomp u(0,0);
 if(l==0){
    for(int j=0; j<Ng;j++){u+=0.5*pi*wg[j]*(ft_cofs[j][0]-0.5*ft_cofs[j][2])*chebU(m,xg[j]);}
 }else{
    for(int j=0; j<Ng;j++){u+=0.25*pi*wg[j]*(ft_cofs[j][l]-ft_cofs[j][l+2])*chebU(m,xg[j]);}
 }
 return u;
}

// k complejo
dcomp pW(int m, int l ,
        Domain & dom,
        dcomp*** ft_cofs,dcomp** wb,int Ng,double* xg,double* wg, dcomp k,int domt,int doms,int dom_op){

 dcomp u(0,0);
 /////////////////////////////////////////W1/////////////////////////////// (el -k**2 incluido en el kernell)
 if(l==0){
    for(int j=0; j<Ng;j++){u+=0.5*pi*wg[j]*(ft_cofs[3][j][0]-0.5*ft_cofs[3][j][2])*chebU(m,xg[j]);}
 }else{
    for(int j=0; j<Ng;j++){u+=0.25*pi*wg[j]*(ft_cofs[3][j][l]-ft_cofs[3][j][l+2])*chebU(m,xg[j]);}
 }
 if(domt==doms){
     if(l==0){
        // for(int j=0; j<Ng;j++){u+=-pow(k,2)*0.25*wg[j]*(log(2)-0.5*chebT(2,xg[j]))*chebU(m,xg[j]);}
        for(int j=0; j<Ng;j++){u+= k*k*0.25*wg[j]*(log(2)-0.5*chebT(2,xg[j]))*chebU(m,xg[j]);}
     }else{
        // for(int j=0; j<Ng;j++){u+=-pow(k,2)*0.25*wg[j]*(chebT(l,xg[j])/l-chebT(l+2,xg[j])/(l+2))*chebU(m,xg[j]);}
        for(int j=0; j<Ng;j++){u+= k*k*0.25*wg[j]*(chebT(l,xg[j])/l-chebT(l+2,xg[j])/(l+2))*chebU(m,xg[j]);}
     }
 }
 //////////////////////////////////////////W2//////////////////////////////

 int Nc = dom.GetNc();

 double jt,jtp,um,ump,aux;
 if(l==0){
    for(int j=0;j<Ng;j++){
        um=chebU(m,xg[j]);
        ump=chebU_p(m,xg[j]);
//         jt=J(xg[j],domt,dom_op);
//         jtp=Jp(xg[j],domt,dom_op);



        jt = dom.Jac(j,0,domt);

        jtp= dom.JacP(j,0,domt);

        u+=0.5*pi*(-l-1)*wg[j]*ump*ft_cofs[4][j][l+1];
        u+=0.5*pi*(l+1)*wg[j]*um*jtp*pow(jt,-2)*ft_cofs[5][j][l+1];
        u+=-wg[j]*ump*0.5*pi/jt*(ft_cofs[6][j][0]-0.5*ft_cofs[6][j][2]);
        u+=wg[j]*um*jtp*0.5*pi*pow(jt,-2)*(ft_cofs[6][j][0]-0.5*ft_cofs[6][j][2]);
    }
 }else{
    for(int j=0;j<Ng;j++){
        um=chebU(m,xg[j]);
        ump=chebU_p(m,xg[j]);
//         jt=J(xg[j],domt,dom_op);
//         jtp=Jp(xg[j],domt,dom_op);

        jt = dom.Jac(j,0,domt);

        jtp= dom.JacP(j,0,domt);

        u+=0.5*pi*(-l-1)*wg[j]*ump*ft_cofs[4][j][l+1];
        u+=0.5*pi*(l+1)*wg[j]*um*jtp*pow(jt,-2)*ft_cofs[5][j][l+1];
        u+=-wg[j]*ump*pi*0.25/jt*(ft_cofs[6][j][l]-ft_cofs[6][j][l+2]);
        u+=wg[j]*um*jtp*0.25*pi*pow(jt,-2)*(ft_cofs[6][j][l]-ft_cofs[6][j][l+2]);
    }
 }

 if(domt==doms){
      if(l==0){
        for(int j=0;j<Ng;j++){
            um=chebU(m,xg[j]);
            ump=chebU_p(m,xg[j]);
//             jt=J(xg[j],domt,dom_op);
//             jtp=Jp(xg[j],domt,dom_op);
            aux=(log(2)-chebT(2,xg[j])/2);

            jt = dom.Jac(j,0,domt);

            jtp= dom.JacP(j,0,domt);

            u+=wg[j]*ump*(-0.5*chebT(1,xg[j])*pow(jt,-2));
            u+=0.5*wg[j]*um*jtp*chebT(1,xg[j])*pow(jt,-3);
            u+=-wg[j]*ump*jtp*pow(jt,-3)*0.25*aux;
            u+=wg[j]*um*(pow(jtp,2)*pow(jt,-4))*0.25*aux;
        }
     }else{
        for(int j=0;j<Ng;j++){
            um=chebU(m,xg[j]);
            ump=chebU_p(m,xg[j]);
//             jt=J(xg[j],domt,dom_op);
//             jtp=Jp(xg[j],domt,dom_op);
            aux=(chebT(l,xg[j])/l-chebT(l+2,xg[j])/(l+2));

            jt = dom.Jac(j,0,domt);

            jtp= dom.JacP(j,0,domt);

            u+=wg[j]*ump*(-0.5*chebT(l+1,xg[j])/pow(jt,2));
            u+=0.5*wg[j]*um*jtp*chebT(l+1,xg[j])/pow(jt,3);
            u+=-wg[j]*ump*jtp*pow(jt,-3)*0.25*aux;
            u+=wg[j]*um*(pow(jtp,2)*pow(jt,-4))*0.25*aux;
        }
     }
 }
 ////////////////////////////////////W3 (signo menos) ////////////////////////////////////


 //integral de simple de (VU_m(t)/J(t))|(t=-1 ,1)*(w U_l/J(s))'
 // como t toma valores en el borde se pude cmabiar la interfaz para que coinsida con la de s
 //en ese caso se regulariza.

//  int bol=reg(domt,doms,dom_op);

 int bol = dom.reg(domt,doms);

 int r1,r2;
 r1=1;
 r2=0;

 double jM= dom.Jac(Nc-1,1,domt); //J(1,domt,dom_op);
 double jm=dom.Jac(0,1,domt); //J(-1,domt,dom_op);

 u-=wb[r1][l+1]*0.5*pi*(double)(l+1)*(double)(-m-1)/jM;
 u-=wb[r2][l+1]*0.5*pi*(double)(l+1)*(double)(m+1)*pow(-1,m)/jm;

 double js,jsp;

 if(l==0){
    u-=(wb[r1+2][0]-0.5*wb[r1+2][2])*0.5*pi*(double)(-m-1)/jM;
    u-=(wb[r2+2][0]-0.5*wb[r2+2][2])*0.5*pi*(double)(m+1)*pow(-1,m)/jm;
 }else{
    u-=0.5*(wb[r1+2][l]-wb[r1+2][l+2])*0.5*pi*(double)(-m-1)/jM;
    u-=0.5*(wb[r2+2][l]-wb[r2+2][l+2])*0.5*pi*(double)(m+1)*pow(-1,m)/jm;
 }
 //regularizacion

 if((bol==0)||(bol==2)||(bol==1)){
     int r=-1;

     js= dom.Jac(0,1,doms);//jsN;

     jsp = dom.JacP(0,1,doms);//jspN;

     if(bol==2){
         r=1;

         js= dom.Jac(Nc-1,1,doms);//jsP;

         jsp =dom.JacP(Nc-1,1,doms);// jspP;
     }
     u-=-(double)(m+1)*0.5*chebT(l+1,r)/(jM*js);
     if(l==0){
        u-=(log(2)-0.5*chebT(2,r))*0.25*(double)(-m-1)*jsp/(jM*pow(js,2));
     }else{
        u-=(chebT(l,r)/(l)-chebT(l+2,r)/(l+2))*0.25*(double)(-m-1)*jsp/(jM*pow(js,2));
     }
 }
 if((bol==0)||(bol==2)||(bol==-1)){
   int r=1;

   js= dom.Jac(Nc-1,1,doms);//jsP;

   jsp =dom.JacP(Nc-1,1,doms);// jspP;

   if(bol==2){
        r=-1;

        js= dom.Jac(0,1,doms);//jsN;

        jsp = dom.JacP(0,1,doms);//jspN;
   }
   u-=(double)(m+1)*pow(-1,m)*0.5*chebT(l+1,r)/(jm*js);
   if(l==0){
        u-=-(log(2)-chebT(2,r)/2)*0.25*(double)(-m-1)*pow(-1,m)*jsp/(jm*pow(js,2));
   }else{
        u-=-(chebT(l,r)/l-chebT(l+2,r)/(l+2))*0.25*(double)(-m-1)*pow(-1,m)*jsp/(jm*pow(js,2));
   }
 }


 //////////////////////////////////////////////////////////////////////////
 return u;
}


double p1(dcomp* ft_cofs, int l){
    if(l==0){
        return 0.5*pi*real(ft_cofs[0]-0.5*ft_cofs[2]);
    }else{
        return 0.25*pi*real(ft_cofs[l]-ft_cofs[l+2]);
    }
}

dcomp pld(int l,dcomp* ft_cofs,int dom_op){
    int sig= pow(-1,l);
    if(dom_op==0){
        sig=1;
    }
    if(l==0){
        return sig*0.5*pi*(ft_cofs[0]-0.5*ft_cofs[2]);
    }else{
        return sig*0.25*pi*(ft_cofs[l]-ft_cofs[l+2]);
    }
}

// k complejo
dcomp pld_d(int l,int Nc,double xc[],double wc[], int doms, int dom_op,dcomp k,int ll, double R){
//     dcomp u(0,0);
//     int sign=-1;
//     if(dom_op==0){
//         sign=1;
//     }
//     for (int j=0;j<Nc;j++){
//         u+=wc[j]*uinc_d(k,xc[j],doms ,dom_op,ll,R)*chebU(l,sign*xc[j]);
//     }
//     return u;
}

// k complejo
dcomp pld_n(int l,int Nc,double xc[],double wc[], int doms, int dom_op,dcomp k, int ll, double R){
//     dcomp u(0,0);
//     int sign=-1;
//     if(dom_op==0){
//         sign=1;
//     }
//     for (int j=0;j<Nc;j++){
//         u+=wc[j]*uinc_n(k,xc[j],doms ,dom_op,ll,R)*chebU(l,sign*xc[j]);
//     }
//     return u;
}



///////////////////otras cosas...//////////////////////////////////////////

double *t_polynomial ( int m, int n, double x[] ){
  int i;
  int j;
  double *v;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    v[i] = 1.0;
  }
  if ( n < 1 )
  {
    return v;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = x[i];
  }
  for ( j = 2; j <= n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      v[i+j*m] = 2.0 * x[i] * v[i+(j-1)*m] - v[i+(j-2)*m];
    }
  }
  return v;
}


double chebT ( int n, double x ){
  int m;
  double *v_vec;
  double value;
  double x_vec[1];

  m = 1;
  x_vec[0] = x;

  v_vec = t_polynomial ( m, n, x_vec );

  value = v_vec[n];

  delete [] v_vec;

  return value;
}

double *u_polynomial ( int m, int n, double x[] ){
  int i;
  int j;
  double *v;

  if ( n < 0 )
  {
    return NULL;
  }

  v = new double[m*(n+1)];

  for ( i = 0; i < m; i++ )
  {
    v[i+0*m] = 1.0;
  }

  if ( n < 1 )
  {
    return v;
  }

  for ( i = 0; i < m; i++ )
  {
    v[i+1*m] = 2.0 * x[i];
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 2; j <= n; j++ )
    {
      v[i+j*m] = 2.0 * x[i] * v[i+(j-1)*m] - v[i+(j-2)*m];
    }
  }
  return v;
}

double chebU ( int n, double x ){
  int m;
  double *v_vec;
  double value;
  double x_vec[1];

  m = 1;
  x_vec[0] = x;

  v_vec = u_polynomial ( m, n, x_vec );

  value = v_vec[n];

  delete [] v_vec;

  return value;
}


double  chebU_p ( int n, double x ){
  return ((n+1)*chebT(n+1,x)-x*chebU( n, x))/(pow(x,2)-1);
}
