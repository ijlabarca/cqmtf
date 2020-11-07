#include "mex.h"
#include "Adom.h"
#include "geo.h"
#include "kers.h"
#include "greenfunction.h"
#include "domain.h"
#include <omp.h>
#include <complex>
#include "IncidentField.h"
using namespace std;
typedef complex<double> dcomp;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //Which mex function
    int fun;
    //fun=1 ==> A
    //fun=2 ==> X
    //fun=3 ==> b
    //fun=4 ==> j
    fun=(int)mxGetScalar(prhs[0]);
    if(fun==1){
        //inpuut variables
        // k complejo
        // double k;
        double kreal, kimag;
        int dom_op,Ng,Nc,Dt,NT,*N;


        //poiters for arrays coming from matlab.
        double *xg,*wg,*pR,*pI,*N_m;
        mxArray *p_m;

        //get scalar variables.
        kreal = mxGetScalar(prhs[1]);
        kimag = mxGetScalar(prhs[2]);
        dcomp k(kreal, kimag);
        N_m= mxGetPr(prhs[3]);            // Number of Chebyschev polynomials per interface
        Nc=(int)mxGetScalar(prhs[4]);     // Chebyschev points
        Ng=(int)mxGetScalar(prhs[5]);     // Gauss-Legendre quad. points
        dom_op=(int)mxGetScalar(prhs[6]); // Domain trial
        Dt=(int)mxGetScalar(prhs[7]);     // Number of interfaces

        double *GaussData, *ChebData;

        GaussData = mxGetPr(prhs[8]);

        ChebData = mxGetPr(prhs[9]);

        Domain dom(Dt,Ng,Nc,
                GaussData,
                ChebData);

        /////////////////////

        xg=mxGetPr(prhs[10]);
        wg=mxGetPr(prhs[11]);

        N=new int[Dt];
        for (int t=0; t <Dt; t++){
            N[t]=(int)N_m[t]; }

        // Matrix
        NT=0;
        for (int t=0; t<Dt; t++){
            NT+=2*N[t]+1;           // 2N+1 degrees of freedom per interface
        }

        dcomp** A=new dcomp*[2*NT];
        for(int j=0; j<2*NT;j++){
            A[j]= new dcomp[2*NT];
        }


        Adom(A,k,N,Nc,Ng,xg,wg,dom_op,Dt,NT,dom);

        p_m=plhs[0] =mxCreateDoubleMatrix(2*NT,2*NT,mxCOMPLEX);
        pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
        for (int j=0; j<2*NT;j++){
            for (int l=0; l<2*NT;l++){
               pR[l*2*NT+j]=real(A[j][l]);
               pI[l*2*NT+j]=imag(A[j][l]);
        }}

       //free memory and end

        for(int j=0; j<2*NT; j++){
          delete A[j];
        }
        delete [] A;
        delete N;

////////////////////////////////////////////////////////////////////////////////


    }else if(fun==2){

        //input variables
        int dom_op_trial,Dt_trial,dom_op_test,Dt_test,NT_trial,NT_test;
        int *N_trial, *N_test;

        //poiters for arrays coming from matlab.
        double *pR,*N_trial_m,*N_test_m;
        mxArray *p_m;

        //get scalar variables.
        N_trial_m=mxGetPr(prhs[1]);
        N_test_m=mxGetPr(prhs[2]);
        dom_op_trial=(int)mxGetScalar(prhs[3]); // Domain trial
        Dt_trial=(int)mxGetScalar(prhs[4]);     // Number of interfaces trial
        dom_op_test=(int)mxGetScalar(prhs[5]);  // Domain test
        Dt_test=(int)mxGetScalar(prhs[6]);      // Number of interfaces test


        int Ng,Nc;

        double *GaussData, *ChebData;

        GaussData = mxGetPr(prhs[7]);

        ChebData = mxGetPr(prhs[8]);

        Ng=(int)mxGetScalar(prhs[9]);

        Nc=(int)mxGetScalar(prhs[10]);

        Domain domTrial(Dt_trial,Ng,Nc,
        GaussData,
        ChebData);


        mxArray *rhs[2];

        rhs[0] = mxDuplicateArray(prhs[11]);

        rhs[1] = mxDuplicateArray(prhs[12]);


        ////////////

        N_trial=new int[Dt_trial];
        for (int t=0; t <Dt_trial; t++){
            N_trial[t]=(int)N_trial_m[t]; }
        N_test=new int[Dt_test];
        for (int t=0; t <Dt_test; t++){
            N_test[t]=(int)N_test_m[t]; }

        //generate matrix.
        NT_trial=0;
        NT_test=0;
        for (int t=0; t<Dt_trial; t++){
            NT_trial+=2*N_trial[t]+1;
        }
        for (int t=0; t<Dt_test; t++){
            NT_test+=2*N_test[t]+1;
        }

        double** X=new double*[2*NT_test];
        for(int j=0; j<2*NT_test;j++){
            X[j]= new double[2*NT_trial];
            for (int l=0 ; l<2*NT_trial; l++){
             X[j][l]=0;
            }
        }

        Xdoms(X,N_trial,N_test,dom_op_trial,Dt_trial,dom_op_test,
                Dt_test,NT_trial,NT_test,domTrial,rhs);

        //output variable
        p_m=plhs[0] =mxCreateDoubleMatrix(2*NT_test,2*NT_trial,mxREAL);
        pR=  mxGetPr(p_m);
        for (int j=0; j<2*NT_test;j++){
            for (int l=0; l<2*NT_trial;l++){
                pR[l*2*NT_test+j]=X[j][l];
        }}

        //free memory and end
        for(int j=0; j < 2*NT_test; j++)
        {
          delete X[j];
        }
        delete [] X;
        delete N_test;
        delete N_trial;

        mxDestroyArray(rhs[0]);
        mxDestroyArray(rhs[1]);



    }else if(fun==3){
       //inpuut variables
       // k complejo
        double R;
        dcomp k0;
        int dom_op,NT,Dt,Dt_0,l;
        double *pR,*pI,*N_m;
        mxArray *p_m;
        int *N;
        //get scalar variables.
        N_m=mxGetPr(prhs[1]);
        dom_op=(int)mxGetScalar(prhs[2]);
        Dt=(int)mxGetScalar(prhs[3]);
        Dt_0=(int)mxGetScalar(prhs[4]);
       ///Domain inicializaiton/////

        int Ng0,Nc0;

        double *GaussData0, *ChebData0;

        GaussData0 = mxGetPr(prhs[5]);

        ChebData0 = mxGetPr(prhs[6]);

        Ng0=(int)mxGetScalar(prhs[7]);
//
        Nc0=(int)mxGetScalar(prhs[8]);

        Domain dom0(Dt_0,Ng0,Nc0,
        GaussData0,
        ChebData0);

        mxArray *rhs[2];

        rhs[0] = mxDuplicateArray(prhs[9]);

        rhs[1] = mxDuplicateArray(prhs[10]);

        mxArray *inc[1];

        inc[0] = mxDuplicateArray(prhs[11]);

        IncidentField incident(dom0,inc);

        N=new int[Dt];
        for (int t=0; t <Dt; t++){
            N[t]=(int)N_m[t]; }


        //generete vector
        NT=0;
        for (int t=0; t<Dt; t++){
            NT+=2*N[t]+1;
        }
        dcomp* b=new dcomp[2*NT];
        for (int j=0; j<2*NT; j++){
            b[j]=dcomp(0,0);
        }

        bdom(b,N,dom_op,Dt,Dt_0,incident,NT,dom0,rhs);

        //output variable
        p_m=plhs[0] =mxCreateDoubleMatrix(2*NT,1,mxCOMPLEX);
        pR=  mxGetPr(p_m);
        pI=  mxGetPi(p_m);
        for (int j=0; j<2*NT;j++){
           pR[j]=real(b[j]);
           pI[j]=imag(b[j]);
        }

        //free memory and end
        delete [] b;
        delete N;

        mxDestroyArray(rhs[0]);
        mxDestroyArray(rhs[1]);



    }

    return;
}
