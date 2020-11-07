% Clase of a interface, it stores an open curve parametrized in [-1,1]
% used as boundarie of domains.  



classdef Interface
    
        
    properties
        m_id
        m_x
        m_y
        m_xp
        m_yp
        m_xpp
        m_ypp
        m_direction
        m_exterior
    end
        
    methods
        
        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = Interface(x,y,dx,dy,ddx,ddy,isExterior,direction,id)
            
             if nargin< 7
               
                isExterior = 1; 
                
            end
            
            if nargin < 8
                
                direction = 1;
                
            end   
            
            if nargin < 9
                
                id = 28;
                
            end   
            
            self.m_x =x; 
            
            self.m_y = y; 
            
            self.m_xp = dx;
            
            self.m_yp = dy;
            
            self.m_xpp = ddx;
            
            self.m_ypp = ddy;
            
            self.m_direction = direction; 
            
            self.m_id = id;
            
            self.m_exterior=isExterior;
        
        end
        
        function [x,y] = geo(self,t)         
            
            x = self.m_x( self.m_direction*t); 
            
            y = self.m_y( self.m_direction*t); 
            
        end 
        
        function [x,y] = geop(self,t)         
            
            x = self.m_direction*self.m_xp( self.m_direction*t); 
            
            y = self.m_direction*self.m_yp( self.m_direction*t); 
            
        end 
        
         function [x,y] = geopp(self,t)         
            
            x = self.m_xpp( self.m_direction*t); 
            
            y = self.m_ypp( self.m_direction*t); 
            
         end 
         
         function u =J(self,t)
             
             [x,y] = self.geop(t); 
             
             u = sqrt(x.*x+y.*y);
             
         end 
         
         function u = Jp(self,t)
             
             [xp,yp] = self.geop(t); 
             
             [xpp,ypp] = self.geopp(t); 
             
             u = (xp.*xpp+yp.*ypp)./sqrt(xp.*xp+yp.*yp);       
             
         end 
        
         function [nx,ny] = normal(self,t)
   
            [xp,yp] = self.geop(t); 
            
            j= sqrt(xp.*xp+yp.*yp);
            
            nx = yp./j; 
             
            ny = -xp./j; 
             
         end 
                  
         function self=SetDirection(self,dir)
             
             self.m_direction =dir;
             
         end 
        
    end 
    
    
end