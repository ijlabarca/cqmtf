#ifndef MY_HEADER_H3
#define MY_HEADER_H3
#include "mex.h"
#include "domain.h"
#include <complex>  

using namespace std;

typedef complex<double> dcomp;

class IncidentField {
    
protected: 
    
    Domain* m_dom;
   
    mxArray *m_incident[1];   
            
    int m_Nc;  
    
    int m_NIntrz;
    
    dcomp* m_DirichletTrace; 
    
    dcomp* m_NeumanTrace; 
    
public:

    IncidentField(Domain &dom,mxArray *incident[1]):
    m_dom(&dom)  
    {
        m_incident[0] = incident[0];
        
        m_Nc = dom.GetNc();
        
        m_NIntrz = dom.GetNItrfz();
        
        m_DirichletTrace =  new dcomp[m_NIntrz*m_Nc]; 
        
        m_NeumanTrace = new dcomp[m_NIntrz*m_Nc];       
        
        Setup(); 
        
    }               
        
    ~IncidentField()  
    {
        delete m_DirichletTrace; 
        
        delete  m_NeumanTrace; 
    }    
    
   void Setup(); 
    
   dcomp DirichletTraze(int n, int intrfz); 
   
   dcomp NeumanTraze(int n, int intrfz); 

   //functions  that call matlab,each time, very slow 
   //and shouldn't be used. 
   dcomp DirichletTraze(double t, int intrfz); 
   
   dcomp NeumanTraze(double s, int intrfz); 
   
        
}; 


#endif