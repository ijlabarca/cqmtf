#include "kers.h"
#include "prods.h"
#include "geo.h"
#include "tchv_cof.h"
#include "in_cells.h"
#include <omp.h>
#include <fftw3.h>
#include "domain.h"
#include "greenfunction.h"
#include "mex.h"

using namespace std;
typedef complex<double> dcomp;



// k complejo
void in_cells(dcomp k,dcomp***** Opd,dcomp**** Wb,int Dt, int dom_op, int Ng, double*xg, int Nc,double* xc,double xp, Domain& dom){
    fftw_plan p;
    int nthread = omp_get_max_threads();
    omp_set_num_threads(nthread);


    //sacar la paralelizacion para trabajr con geometrias de matlab?

    #pragma omp parallel
    {

        FreeSpace gf(k,&dom);

        #pragma omp for collapse(4)
        for(int t=0; t<Dt;t++){
            for(int s=0; s<Dt;s++){
                for(int j=0; j<Ng;j++){
                    for(int l=0;l<Nc;l++){
                        Opd[t][s][0][j][l]=GB(xg[j], xc[l],
                                j,l,0,1,t,s, dom_op,gf);
                        Opd[t][s][1][j][l]=KB(j,l,0,1,t,s, dom_op,gf);
                        Opd[t][s][2][j][l]=KAB(j,l,0,1,t,s, dom_op,gf);
                        Opd[t][s][3][j][l]=WB(xg[j], xc[l],
                                j,l,0,1,t,s, dom_op,gf,dom);
                        Opd[t][s][4][j][l]=WA1(xg[j], xc[l],
                                j,l,0,1,t,s, dom_op,gf,dom);
                        Opd[t][s][5][j][l]=WA2(xg[j], xc[l],
                                j,l,0,1,t,s, dom_op,gf,dom);
                        Opd[t][s][6][j][l]=WA3(xg[j], xc[l],
                                j,l,0,1,t,s, dom_op,gf,dom);
                        if(j==0){
                            Wb[t][s][0][l]= WB1(-1,xg[Ng-1],xc[l],
                                    Ng-1,l,0,1,t,s,dom_op,gf,dom);

                            Wb[t][s][1][l]= WB1(1,xg[0],xc[l],
                                    0,l,0,1,t,s,dom_op,gf,dom);

                            Wb[t][s][2][l]= WB2(-1,xg[Ng-1],xc[l],
                                    Ng-1,l,0,1,t,s,dom_op,gf,dom);

                            Wb[t][s][3][l]= WB2(1,xg[0],xc[l],
                                    0,l,0,1,t,s,dom_op,gf,dom);
                        }
                    }
                }
            }
        }
    }
       p=plan_0(Nc, Opd[0][0][0][0]  );
       #pragma omp parallel for collapse(2)
       for(int t=0; t<Dt;t++){
           for (int s=0; s< Dt;s++){
                Tchv_cof(Nc,Wb[t][s][0],p);
                Tchv_cof(Nc,Wb[t][s][1],p);
                Tchv_cof(Nc,Wb[t][s][2],p);
                Tchv_cof(Nc,Wb[t][s][3],p);
            }
       }

     #pragma omp parallel for collapse(3)
       for(int t=0; t<Dt;t++){
         for(int s=0; s<Dt;s++){
            for(int j=0; j<Ng;j++){
                Tchv_cof(Nc,Opd[t][s][0][j],p);
                Tchv_cof(Nc,Opd[t][s][1][j],p);
                Tchv_cof(Nc,Opd[t][s][2][j],p);
                Tchv_cof(Nc,Opd[t][s][3][j],p);
                Tchv_cof(Nc,Opd[t][s][4][j],p);
                Tchv_cof(Nc,Opd[t][s][5][j],p);
                Tchv_cof(Nc,Opd[t][s][6][j],p);

    }}}
    fftw_destroy_plan(p);
    fftw_cleanup();
    return;

    //ESTA FUNCION SE PODRIA MEJORAR, como se necesitan guardar arrays
    // de tamaño Nc, bastan los primeros 2N+1 coeficientes.
}


///////////////////////////////////////////////////////////////////////////
void in_cells_d(dcomp*** Opd,int Dt, int dom_op_trial,int dom_op_test,int Nc,double* xc,int* N,Domain &domTrial){
    fftw_plan p;
    int nthread = omp_get_max_threads();
    omp_set_num_threads(nthread);
    int sign=-1;
    if (dom_op_trial==dom_op_test){
        sign=1; }
    double jt;
    //#pragma omp parallel for collapse(3)
    for(int t=0; t<Dt;t++){
        for(int m=-N[t]; m<N[t]+1;m++){
            for(int l=0;l<Nc;l++){

                    if (sign ==1)
                    {jt = domTrial.Jac(l,1,t);}
                    else
                    {jt = domTrial.Jac(Nc-1-l,1,t);}

                    Opd[t][m+N[t]][l]=chebU(m+N[t],(double)sign*xc[l])/jt;
            }
        }
    }

       p=plan_0(Nc, Opd[0][0]);
       //#pragma omp parallel for collapse(2)
       for(int t=0; t<Dt;t++){
           for (int m=-N[t]; m< N[t]+1;m++){
                Tchv_cof(Nc,Opd[t][m+N[t]],p);
            }
       }

    fftw_destroy_plan(p);
    fftw_cleanup();
    return;
    //saca muchas mas transformadas que las que senecesitan
}
//////////////////////////////////////////////////////////////////////////

void in_cells_ld(dcomp*** Opd,int Dt, int Nc,double* xc,IncidentField& inc,Domain&dom0){
    fftw_plan p;
    int nthread = omp_get_max_threads();
    omp_set_num_threads(nthread);



    #pragma omp parallel for collapse(2)
    for(int t=0; t<Dt;t++){
        for(int l=0;l<Nc;l++){
                    Opd[0][t][l]=inc.DirichletTraze(l,t);
//                            %uinc_d(k,l,1,t ,0,ll,R,dom0);
//                     Opd[1][t][l]=0.0;
                    Opd[1][t][l]=inc.NeumanTraze(l,t);

                            //uinc_n(k,l,1,t ,0,ll,R,dom0);
        }
    }


   p=plan_0(Nc, Opd[0][0]);
   #pragma omp parallel for collapse(1)
   for(int t=0; t<Dt;t++){
            Tchv_cof(Nc,Opd[0][t],p);
            Tchv_cof(Nc,Opd[1][t],p);
   }

   fftw_destroy_plan(p);
   fftw_cleanup();
   return;
   //saca muchas mas transformadas que las que senecesitan
}

// k complejo
void in_cells_ldPlaneWave(dcomp*** Opd,int Dt, int Nc,double* xc,
        dcomp k, double alpha,Domain&dom0){
    fftw_plan p;
    int nthread = omp_get_max_threads();
    omp_set_num_threads(nthread);



    #pragma omp parallel for collapse(2)
    for(int t=0; t<Dt;t++){
        for(int l=0;l<Nc;l++){
                    Opd[0][t][l]=planeWave(k,l,1,t ,alpha,dom0);

                    Opd[1][t][l]=planeWaveN(k,l,1,t ,alpha,dom0);
        }
    }


   p=plan_0(Nc, Opd[0][0]);
   #pragma omp parallel for collapse(1)
   for(int t=0; t<Dt;t++){
            Tchv_cof(Nc,Opd[0][t],p);
            Tchv_cof(Nc,Opd[1][t],p);
   }

   fftw_destroy_plan(p);
   fftw_cleanup();
   return;
   //saca muchas mas transformadas que las que senecesitan
}


// k complejo
void in_cells0(dcomp*** uin,int Dt,dcomp k, double alpha, int Nc,double* xc){
    fftw_plan p;
    fftw_destroy_plan(p);
    fftw_cleanup();
}
