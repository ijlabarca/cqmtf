#include <omp.h>
#include "prods.h"
#include "in_cells.h"
#include "tchv_cof.h"
#include "mex.h"
#include "domain.h"
#include "Adom.h"
using namespace std;
typedef complex<double> dcomp;
const double pi= 4.0*atan(1.0);

// k complejo
void Adom(dcomp** A, dcomp k, int* N,int Nc, int Ng, double*xg, double* wg,int dom_op, int Dt, int NT, Domain &dom){


    ///////////////guardamos los coeficientes de fourier para este domonio

    dcomp***** Opd =new dcomp****[Dt];
    for(int t=0;t< Dt; t++){
        Opd[t]=new dcomp***[Dt];
        for(int s=0; s<Dt; s++){
            Opd[t][s]=new dcomp**[8];
            for(int o=0;o<8;o++){
                 Opd[t][s][o]=new dcomp*[Ng];
                 for(int j=0; j<Ng;j++){
                     Opd[t][s][o][j]=new dcomp[Nc];
    }}}}


    dcomp**** Wb= new dcomp***[Dt];
    for(int t=0;t<Dt;t++){
        Wb[t]=new dcomp**[Dt];
        for(int s=0;s<Dt;s++){
            Wb[t][s]=new dcomp*[4];
             for(int o=0;o<4;o++){
                Wb[t][s][o]=new dcomp [Nc];
    }}}

    double* xc;
    xc = chebpts(Nc);
   double xp= 0.99999999999;


   in_cells(k,Opd,Wb,Dt,dom_op,Ng, xg, Nc,xc,xp,dom);

    int ns=0;
    int nt=0;
    for(int t=0;t<Dt;t++){
        if(t>0) nt+=2*N[t-1]+1;
        ns=0;
        for(int s=0;s<Dt;s++){
            if(s>0) ns+=2*N[s-1]+1;
            #pragma omp parallel for collapse(2)
            for(int l=0;l<2*N[s]+1;l++){
                for(int m=-N[t];m<N[t]+1;m++){
                     A[l+ns][m+N[t]+nt]= -2.0*pK(m+N[t],l ,Opd[t][s][1],Ng,xg,wg,k,t,s);
                     A[l+ns][m+N[t]+NT+nt]=2.0*pV(m+N[t],l ,Opd[t][s][0],Ng,xg,wg,k,t,s);
                     A[l+NT+ns][m+N[t]+nt]=2.0* pW(m+N[t],l,dom,
                             Opd[t][s],Wb[t][s],Ng,xg,wg,k,t,s,dom_op);
                     A[l+NT+ns][m+N[t]+NT+nt]=2.0*pK(m+N[t],l ,Opd[t][s][2],Ng,xg,wg,k,t,s);
                 }}}}


     delete [] xc;
     xc=0;



    for(int t=0;t< Dt; t++){
        for(int s=0; s<Dt; s++){
            for(int o=0;o<8;o++){
                 for(int j=0; j<Ng;j++){
                     delete [] Opd[t][s][o][j];}
                 delete []  Opd[t][s][o];}
            delete [] Opd[t][s];}
        delete [] Opd[t];}
    delete [] Opd;
     Opd= 0;






    for(int t=0;t<Dt;t++){
        for(int s=0;s<Dt;s++){
             for(int o=0;o<4;o++){
                delete [] Wb[t][s][o];}
            delete [] Wb[t][s]; }
        delete [] Wb[t];}
    delete [] Wb;
    Wb=0;
    return;


}


void Xdoms(double** X,  int* N_trial,int* N_test,int dom_op_trial,
        int Dt_trial,int dom_op_test, int Dt_test,
        int NT_trial, int NT_test,
        Domain &domTrial,mxArray *doms[2]){
    //trial_test un vector de dimension Dt_trial, sus componentes indican con que subdominios de dom_op_test se comunican
    // los dominios de dom_op_trial, ej trial_test[j]= k indica que el subdominio j se comunica con el subdominio k, si k=-1 no se comunica con ninguno.
    // X inicialisarla en 0.



    int* tt=trial_test(Dt_trial, doms);

    int Nc=domTrial.GetNc();
    int s;
    dcomp*** Opd =new dcomp**[Dt_trial];
    for(int t=0;t< Dt_trial; t++){
        Opd[t]=new dcomp*[2*N_trial[t]+1];
         for(int m=-N_trial[t]; m<N_trial[t]+1;m++){
           Opd[t][m+N_trial[t]]=new dcomp[Nc];
    }}

    double* xc;
    xc = chebpts(Nc);
    in_cells_d(Opd,Dt_trial,dom_op_trial,dom_op_test,Nc,xc,N_trial,domTrial);
    int sign=-1;
    if (dom_op_trial==dom_op_test) {
        sign=1;
        tt=trial(Dt_trial);
    }
    int nt=0;
    int ns;

    for (int t=0;t<Dt_trial;t++){
        if(t>0) nt+=2*N_trial[t-1]+1;
        if(tt[t]>-1){
            s=tt[t];
            ns=0;
            if(s>0){
                for (int as=0;as<s;as++) ns+=2*N_test[as]+1;
            }
            for(int m=-N_trial[t];m<N_trial[t]+1;m++){
                for(int l=0;l<2*N_test[s]+1;l++){
                    X[l+ns][m+N_trial[t]+nt]= (double)sign*p1(Opd[t][m+N_trial[t]],l);
                    X[l+NT_test+ns][m+N_trial[t]+NT_trial+nt]=(double)sign*X[l+ns][m+N_trial[t]+nt];
                }}
        }
    }
    delete tt;
        for(int t=0;t< Dt_trial; t++){
            for(int m=0; m<2*N_trial[t]+1;m++){
                delete []Opd[t][m];
            }
            delete [] Opd[t];
        }
    delete [] Opd;

//     delete aux;




}

void bdom(dcomp* b, int* N,int dom_op, int Dt,int Dt_0,IncidentField& inc,
        int NT,Domain &dom0,mxArray *doms[2]){
 // b inicializado en 0.
     int Nc=dom0.GetNc();
     int* tt;
     int sign_dom=1;
     if(dom_op==0){
         sign_dom=-1;
         tt=trial(Dt);
     }else{
         tt=trial_test(Dt, doms);
     }
     double* xc;
     xc = chebpts(Nc);
     dcomp*** Opd =new dcomp**[2];
     Opd[0]= new dcomp*[Dt_0];
     Opd[1]= new dcomp*[Dt_0];
     for(int t=0;t< Dt_0; t++){
            Opd[0][t]=new dcomp[Nc];
            Opd[1][t]=new dcomp[Nc];
      }

     in_cells_ld(Opd,Dt_0,Nc,xc,inc,dom0);
     int ns=0;
     for (int s=0; s<Dt;s++){
         if( s>0) ns+=2*N[s-1]+1;
         if(tt[s]>-1){
            for(int l=0;l<2*N[s]+1;l++){
                b[l+ns]=(double)sign_dom*pld(l,Opd[0][tt[s]],dom_op);
                b[l+ns+NT]=-pld(l,Opd[1][tt[s]],dom_op);
            }
         }
     }
      delete tt;
          for(int t=0;t< Dt_0; t++){
            delete [] Opd[0][t];
            delete [] Opd[1][t];
      }
      delete [] Opd[0];
      delete [] Opd[1];
      delete [] Opd;

      return;
}

// k complejo
void bdomPlaneWave(dcomp* b, int* N,int dom_op, int Dt,int Dt_0,
        dcomp k,double alpha,int NT,Domain &dom0,mxArray *doms[2]){
 // b inicializado en 0.
     int Nc=dom0.GetNc();
     int* tt;
     int sign_dom=1;
     if(dom_op==0){
         sign_dom=-1;
         tt=trial(Dt);
     }else{
         tt=trial_test(Dt, doms);
     }
     double* xc;
     xc = chebpts(Nc);
     dcomp*** Opd =new dcomp**[2];
     Opd[0]= new dcomp*[Dt_0];
     Opd[1]= new dcomp*[Dt_0];
     for(int t=0;t< Dt_0; t++){
            Opd[0][t]=new dcomp[Nc];
            Opd[1][t]=new dcomp[Nc];
      }

     in_cells_ldPlaneWave(Opd,Dt_0,Nc,xc,k,alpha,dom0);
     int ns=0;
     for (int s=0; s<Dt;s++){
         if( s>0) ns+=2*N[s-1]+1;
         if(tt[s]>-1){
            for(int l=0;l<2*N[s]+1;l++){
                b[l+ns]=(double)sign_dom*pld(l,Opd[0][tt[s]],dom_op);
                b[l+ns+NT]=-pld(l,Opd[1][tt[s]],dom_op);
            }
         }
     }
      delete tt;
          for(int t=0;t< Dt_0; t++){
            delete [] Opd[0][t];
            delete [] Opd[1][t];
      }
      delete [] Opd[0];
      delete [] Opd[1];
      delete [] Opd;

      return;
}

int* trial_test(int Dt_trial,mxArray *doms[2])
{
    int* tt = new int[Dt_trial];

    double *aux;

    mxArray *lhs[1];

    mexCallMATLAB(1, lhs, 2, doms, "TrialTest");

    aux = mxGetPr(lhs[0]);

    for(int ii(0);ii<Dt_trial;++ii)
    {
        tt[ii] = (int)aux[ii];

    }

    mxDestroyArray(lhs[0]);

    return tt;

}

int* trial(int Dt){
 int *tt =new int[Dt];
 for(int j=0; j< Dt; j++){
     tt[j]=j;
 }
 return tt;
}
