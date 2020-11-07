% Clase of a interface, it stores an open curve parametrized in [-1,1]
% used as boundarie of domains.



classdef Domain


    properties
        m_id
        m_NumInfefaces
        m_Interfaces
        m_KwaveNumber
    end

    methods

        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = Domain(id,kreal, kimag,Interfaces)

            if nargin< 4

                Interfaces = [];

            end

            self.m_id = id;

            self.m_KwaveNumber = kreal + 1i * kimag;  

            self.m_Interfaces = Interfaces;

            self.m_NumInfefaces = length(Interfaces);


        end


        function self = AddInterface(self,Interface)

            self.m_NumInfefaces = 1+self.m_NumInfefaces;

            self.m_Interfaces{self.m_NumInfefaces} = Interface;

        end
% 
%         function self = set.m_KwaveNumber(self, ks)
%            
%             self.m_KwaveNumber = ks(1) + 1i*ks(2);
%             
%             
%         end
        
        
        function [x,y] = geo(self,t,numInterface)
%             disp(numInterface);

            [x,y] = self.m_Interfaces{numInterface}.geo(t);

        end

        function [x,y] = geop(self,t,numInterface)

            [x,y] = self.m_Interfaces{numInterface}.geop(t);

        end

         function [x,y] = geopp(self,t,numInterface)

            [x,y] = self.m_Interfaces{numInterface}.geopp(t);

         end

         function u =J(self,t,numInterface)

           u = self.m_Interfaces{numInterface}.J(t);

         end

         function u = Jp(self,t,numInterface)

            u = self.m_Interfaces{numInterface}.Jp(t);

         end

         function [nx,ny] = normal(self,t,numInterface)

            [nx,ny] = self.m_Interfaces{numInterface}.normal(t);

         end

         function Display(self,clr)

              if nargin< 3

                clr  = 'blue';

             end

             N = self.m_NumInfefaces;

             ts= linspace(-1,1,100);

             XY = zeros(100*N,2);

             for intr = 1:N

                 for t = 1:100

                    [x,y] = self.geo(ts(t),intr);

                    XY(t+100*(intr-1),:)=[x,y];

                 end
             end

             plot(XY(:,1),XY(:,2),'color',clr);


         end

         function u = GetKwave(self)

             u = self.m_KwaveNumber;

         end


%          function u=GetNumIntefs(self)
%
%              u = self.m_NumInfefaces;
%
%          end
%
%        function u=GetNumIntefs(self)
%
%              u = self.m_NumInfefaces;
%
%          end

    end


end
