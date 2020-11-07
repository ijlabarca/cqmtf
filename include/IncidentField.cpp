#include "IncidentField.h"

void IncidentField::Setup()
{
    double *auxR,*auxdxR,*auxdyR; 
    
    double *auxImg,*auxdxImg,*auxdyImg; 
       
    double *x ,*y;
    
    double *nx = new double[m_Nc*m_NIntrz];
    
    double *ny = new double[m_Nc*m_NIntrz];
    
    mxArray *rhs[2]; 
    
    rhs[0] = m_incident[0];
    
    rhs[1] = mxCreateDoubleMatrix(m_Nc*m_NIntrz,1,mxCOMPLEX);
    
    x = mxGetPr(rhs[1]);
    
    y = mxGetPi(rhs[1]);
    
    double *zeros = new double[m_Nc*m_NIntrz];
    
    //could be better if we can get a pointer from the 
    //domain but is not the same structure of data.  
    
    for(int jj(0); jj<m_NIntrz; ++jj)
    {
        for (int ll(0) ; ll< m_Nc; ++ll)
        {
            m_dom->Geo(x+jj*m_Nc+ll,y+jj*m_Nc+ll,ll,1,jj);
            
            m_dom->Normal(nx+jj*m_Nc+ll,ny+jj*m_Nc+ll,ll,1,jj);
            
            zeros[jj*m_Nc+ll] =0; 
        }
    }
    
    mxArray *lhs[1];                
    
    mexCallMATLAB(1, lhs, 2, rhs, "evaluate");
      
    auxR = mxGetPr(lhs[0]);
    
    auxImg = mxGetPi(lhs[0]);
    
    if(!auxImg)
    {
        auxImg = zeros; 
    }
    
    mxArray *lhs2[2]; 
    
    mexCallMATLAB(2, lhs2, 2, rhs, "evaluateGradient");    
    
    auxdxR = mxGetPr(lhs2[0]);
    
    auxdxImg = mxGetPi(lhs2[0]);
    
    if(!auxdxImg)
    {
        auxdxImg = zeros; 
    }
    
    auxdyR = mxGetPr(lhs2[1]);
    
    auxdyImg = mxGetPi(lhs2[1]);
    
    if(!auxdyImg)
    {
        auxdyImg = zeros; 
    }        
    
    for(int jj(0); jj<m_NIntrz; ++jj)
    {
        for (int ll(0) ; ll< m_Nc; ++ll)
        {        
            
            m_DirichletTrace[jj*m_Nc+ll] = dcomp(auxR[jj*m_Nc+ll],
                    auxImg[jj*m_Nc+ll]); 
            
            m_NeumanTrace[jj*m_Nc+ll] = dcomp(auxdxR[jj*m_Nc+ll]*
                    nx[jj*m_Nc+ll]+auxdyR[jj*m_Nc+ll]*ny[jj*m_Nc+ll],
                    auxdxImg[jj*m_Nc+ll]*nx[jj*m_Nc+ll]+
                    auxdyImg[jj*m_Nc+ll]*ny[jj*m_Nc+ll]);
        }
    }
    
    //mxDestroyArray(rhs[0]);
    
    mxDestroyArray(rhs[1]);
    
    mxDestroyArray(lhs[0]);
    
    mxDestroyArray(lhs2[0]);
    
    mxDestroyArray(lhs2[1]);
    
    delete [] nx;
    delete [] ny;
    delete [] zeros; 
    
}

dcomp IncidentField::DirichletTraze(int n, int intrfz)
{
    return m_DirichletTrace[intrfz*m_Nc+n]; 
}
   
dcomp IncidentField::NeumanTraze(int n, int intrfz)
{
    return m_NeumanTrace[intrfz*m_Nc+n]; 
    
}

//functions  that call matlab,each time, very slow 
//and shouldn't be used. 
dcomp IncidentField::DirichletTraze(double t, int intrfz) 
{
    return 0; 
    
}

dcomp IncidentField::NeumanTraze(double s, int intrfz)
{
    
    return 0; 
}