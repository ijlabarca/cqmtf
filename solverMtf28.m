classdef solverMtf28 < solver

    properties

        m_DOMS

        m_DomainArray

        m_AllGeoData

        m_kext_real

        m_kext_imag

        m_kext

        m_kwns_real

        m_kwns_imag

        m_kwns

        m_wavespeed

        m_Ng

        m_Nc

        m_N

        m_A

        m_Ai

        m_M

        m_b

        m_tu

        m_IsSetup

        m_setupB


    end

    methods


        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = solverMtf28(c_0, kwave_real, kwave_imag,incidentField,InterDomains)

            self = self@solver(kwave_real, kwave_imag,incidentField);

            self.m_wavespeed = c_0;

            self.m_DomainArray = InterDomains;

            self.m_kext_real = kwave_real;

            self.m_kext_imag = kwave_imag;

            self.m_kext = kwave_real + 1i*kwave_imag;

            self.m_IsSetup = 0;

            self.m_setupB = 0;

            if (~isempty(InterDomains))

                self.m_DOMS = length(InterDomains);

                self = ConstructDomain0(self);

            end

        end

        %===============================================================
        % these methods must be provided
        %===============================================================

        %-----------------------------------------
        % setup
        %-----------------------------------------

        % These methods set up your solver, eg assembles discretisation
        % matrices etc



        function setup(self, Nmax)

           if(self.m_IsSetup == 0)


            mex -v mex_funs.cpp Adom.cpp kers.cpp prods.cpp geo.cpp...
                in_cells.cpp tchv_cof.cpp greenfunction.cpp IncidentField.cpp...
                domain.cpp ...
                CXXFLAGS="\$CXXFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"...
                -lcomplex_bessel -lgfortran -lfftw3 -lm ;



           end



           self.m_Nc=zeros(1,self.m_DOMS+1)+1000; % + 28*Nmax/2;%+ceil(2.8*maxK);

           self.m_Ng=zeros(1,self.m_DOMS+1)+1000; %+ 28*Nmax/2;%+ceil(2.8*maxK);

           %store all geometry.
           self.m_AllGeoData = cell(self.m_DOMS+1,1);

           domt = zeros(self.m_DOMS+1,1);

           self.m_N = cell(self.m_DOMS+1,1);

           for i=1:(self.m_DOMS+1)

               domt(i) = self.m_DomainArray{i}.m_NumInfefaces;

               % Número de grados de libertad
                self.m_N{i} = zeros(domt(i),1)+Nmax;



               xc = cos( (self.m_Nc(i)-1-(0:self.m_Nc(i)-1))*pi...
                   /(self.m_Nc(i)-1));

               [xg,wg]=lgwt(self.m_Ng(i),-1,1);

               self.m_AllGeoData{i} = SetupDoms( self.m_DomainArray,...
                   i-1,xg,xc );


                   %DGauss=[XG;YG;NxG;NyG;JG;JpG];

                   %DCheb =[XC;YC;NxC;NyC;JC;JpC];

                   %AllGeoData{1} = DGauss;

                   %AllGeoData{2} = DCheb;

           end

        dim=0;

        for d=0:self.m_DOMS

            for t=1:domt(d+1)

                dim=dim+2*(2*self.m_N{d+1}(t)+1);

            end
        end

        self.m_A = zeros(dim,dim);

        self.m_Ai = zeros(dim,dim);

        self.m_M = zeros(dim,dim);

        if self.m_setupB == 1
        %% Compute matrix.
        % d_test: number of domain test that is computed
        % d_trial: number of domain trial that is computed
        %
        % for given d_test, d_trial it has to be determined what position of
        % the matrix is begin computed, for that the variables (s_test, s_trial
        % NT_trial, N_test) are used
        %
        %
        %s_test/_trial: gives the start position of the test/trial block.
        %NT_test/trial: gives the number of degrees of fredom for one of the
        %traces (Dirichlet and Neumann use the same number) for a given domain.
        s_test=1;

        for d_test=0:self.m_DOMS

            NT_test=0;

            if(d_test>0)
                %compute the start position(on the matrix) of the test domain
                for t=1:domt(d_test)

                    s_test=s_test+2*(2*self.m_N{d_test}(t)+1);

                end
            end

            %compute the degrees of fredom for the test domain

            for t=1:domt(d_test+1)

                NT_test=NT_test+(2*self.m_N{d_test+1}(t)+1);

            end

            s_trial=1;

            for d_trial=0:self.m_DOMS

                NT_trial=0;
                if(d_trial>0)
                 %compute the start position(on the matrix) of the trial domain
                    for t=1:domt(d_trial)

                        s_trial=s_trial+2*(2*self.m_N{d_trial}(t)+1);

                    end
                end
                %compute the degrees of fredom for the trial domain
                for t=1:domt(d_trial+1)

                    NT_trial=NT_trial+(2*self.m_N{d_trial+1}(t)+1);

                end
               %computation of the diagonal blocks
               if(d_trial==d_test)
                    %gaus legendre points for the given domain
                    [xg,wg]=lgwt(self.m_Ng(d_trial+1),-1,1);



                    %Integral operators block

                        self.m_A( s_test:(s_test-1)+2*NT_test,...
                            s_trial:(s_trial-1)+2*NT_trial)=...
                            ...
                        mex_funs(1,...
                        self.m_kwns_real(d_trial+1),...
                        self.m_kwns_imag(d_trial+1),...
                        self.m_N{d_trial+1},...
                        self.m_Nc(d_trial+1),...
                        self.m_Ng(d_trial+1),...
                        d_trial,...
                        domt(d_trial+1), ...
                        self.m_AllGeoData{d_trial+1}{1},...   % Gauss Data
                        self.m_AllGeoData{d_trial+1}{2},...   % Cheb. Data
                        xg,wg);                               % Gauss-Legendre Quadrature

                    %Duplication of the integral operators block
                    self.m_Ai( s_test:(s_test-1)+2*NT_test,...
                        s_trial:(s_trial-1)+2*NT_trial)=...
                        ...
                         self.m_A( s_test:(s_test-1)+2*NT_test, ...
                         s_trial:(s_trial-1)+2*NT_trial);



                     % This only has to be computed once, not for all wavenumbers
                     % Every call to mex_funs is expensive
                     self.m_M( s_test:(s_test-1)+2*NT_test,...
                         s_trial:(s_trial-1)+2*NT_trial) =...
                         ...
                         mex_funs(2,self.m_N{d_trial+1},self.m_N{d_test+1},...
                         d_trial,domt(d_trial+1),d_test,domt(d_test+1),...
                         self.m_AllGeoData{d_trial+1}{1},...
                         self.m_AllGeoData{d_trial+1}{2},...
                         self.m_Ng(d_trial+1),self.m_Nc(d_trial+1),...
                         self.m_DomainArray{d_trial+1},...
                         self.m_DomainArray{d_trial+1});

               else

               % This only has to be computed once, not for all wavenumbers
               % Every call to mex_funs is expensive
                   self.m_A( s_test:(s_test-1)+2*NT_test,...
                       s_trial:(s_trial-1)+2*NT_trial)=...
                       ...
                       mex_funs(2,self.m_N{d_trial+1},self.m_N{d_test+1},...
                       d_trial,domt(d_trial+1),d_test,domt(d_test+1),...
                        self.m_AllGeoData{d_trial+1}{1},...
                        self.m_AllGeoData{d_trial+1}{2},...
                        self.m_Ng(d_trial+1),self.m_Nc(d_trial+1),...
                        self.m_DomainArray{d_trial+1},...
                        self.m_DomainArray{d_test+1});
               end
            end


        end

        self.m_IsSetup = 1;

        end

        end


        %-----------------------------------------
        % setupb
        %-----------------------------------------

        % This method constructs the matrix rhs. incidentField corresponds
        % to

        function setupb(self, Nmax, M)
%            Nrhs = length(self.incidentField);

            Nrhs = M+1;

            self.setup(Nmax);


            self.m_setupB = 1;

           domt = zeros(self.m_DOMS+1);

           for i=1:(self.m_DOMS+1)

               domt(i) = self.m_DomainArray{i}.m_NumInfefaces;

           end

           dim=0;

           for d=0:self.m_DOMS

               for t=1:domt(d+1)

                    dim=dim+2*(2*self.m_N{d+1}(t)+1);
%                       dim=dim+2*(2*Nmax+1);

               end

           end

           self.m_b = zeros(dim,Nrhs);


           s_test =1;

           for ii =1:Nrhs
%
%                disp(ii);

               s_test =1;

              for d_test=0:self.m_DOMS

                   NT_test=0;

                   if(d_test>0)
                        %compute the start position(on the matrix) of the test domain
                       for t=1:domt(d_test)

                           s_test=s_test+2*(2*self.m_N{d_test}(t)+1);
%                            s_test=s_test+2*(2*Nmax+1);

                       end
                   end
                    %compute the degrees of fredom for the test domain
                   for t=1:domt(d_test+1)

                       NT_test=NT_test+(2*self.m_N{d_test+1}(t)+1);
%                          NT_test=NT_test+(2*Nmax+1);

                   end



                  self.incidentField{1}.m = ii;
                  self.m_b(s_test:(s_test-1)+2*NT_test,ii) = ...
                  mex_funs(3,self.m_N{d_test+1},d_test,domt(d_test+1),...
                           domt(1),self.m_AllGeoData{1}{1},...
                           self.m_AllGeoData{1}{2},self.m_Ng(1),self.m_Nc(1),...
                           self.m_DomainArray{d_test+1},self.m_DomainArray{1},...
                           self.incidentField{1});

              end
          end


        end






        function setupb_freq(self, Nmax)
           Nrhs = length(self.incidentField);



            self.setup(Nmax);

            self.m_setupB = 1;

           domt = zeros(self.m_DOMS+1);

           for i=1:(self.m_DOMS+1)

               domt(i) = self.m_DomainArray{i}.m_NumInfefaces;

           end

           dim=0;

           for d=0:self.m_DOMS

               for t=1:domt(d+1)

                    dim=dim+2*(2*self.m_N{d+1}(t)+1);
%                       dim=dim+2*(2*Nmax+1);

               end

           end

           self.m_b = zeros(dim,Nrhs);


           s_test =1;

           for ii =1:Nrhs
%
%                disp(ii);

               s_test =1;

              for d_test=0:self.m_DOMS

                   NT_test=0;

                   if(d_test>0)
                        %compute the start position(on the matrix) of the test domain
                       for t=1:domt(d_test)

                           s_test=s_test+2*(2*self.m_N{d_test}(t)+1);
%                            s_test=s_test+2*(2*Nmax+1);

                       end
                   end
                    %compute the degrees of fredom for the test domain
                   for t=1:domt(d_test+1)

                       NT_test=NT_test+(2*self.m_N{d_test+1}(t)+1);
%                          NT_test=NT_test+(2*Nmax+1);

                   end


                  self.m_b(s_test:(s_test-1)+2*NT_test,ii) = ...
                  mex_funs(3,self.m_N{d_test+1},d_test,domt(d_test+1),...
                           domt(1),self.m_AllGeoData{1}{1},...
                           self.m_AllGeoData{1}{2},self.m_Ng(1),self.m_Nc(1),...
                           self.m_DomainArray{d_test+1},self.m_DomainArray{1},...
                           self.incidentField{ii});

              end
          end


        end


        %-----------------------------------------
        % solve
        %-----------------------------------------

        % This method solves the scattering problem for every right hand
        % side specified in the self.incidentField cell array.

        function solve(self)

	    Nrhs = length(self.incidentField);

	    self.setup();

           domt = zeros(self.m_DOMS+1);

           for i=1:(self.m_DOMS+1)

               domt(i) = self.m_DomainArray{i}.m_NumInfefaces;

           end

           dim=0;

           for d=0:self.m_DOMS

               for t=1:domt(d+1)

                    dim=dim+2*(2*self.m_N{d+1}(t)+1);

               end

           end

           self.m_b = zeros(dim,Nrhs);

           self.m_tu = zeros(dim,Nrhs);

           s_test =1;

          for ii =1:Nrhs

               s_test =1;

              for d_test=0:self.m_DOMS

                   NT_test=0;

                   if(d_test>0)
                        %compute the start position(on the matrix) of the test domain
                       for t=1:domt(d_test)

                           s_test=s_test+2*(2*self.m_N{d_test}(t)+1);

                       end
                   end
                    %compute the degrees of fredom for the test domain
                   for t=1:domt(d_test+1)

                       NT_test=NT_test+(2*self.m_N{d_test+1}(t)+1);

                   end

                  self.m_b(s_test:(s_test-1)+2*NT_test,ii) = ...
                  mex_funs(3,self.m_N{d_test+1},d_test,domt(d_test+1),...
                           domt(1),self.m_AllGeoData{1}{1},...
                           self.m_AllGeoData{1}{2},self.m_Ng(1),self.m_Nc(1),...
                           self.m_DomainArray{d_test+1},self.m_DomainArray{1},...
                           self.incidentField{ii});

              end
          end

	     self.m_tu = self.m_A\self.m_b;

           %norm(self.m_Ai*self.m_tuplaneWave-self.m_M*...
           %    self.m_tuplaneWave,Inf)




        end


        %===============================================================
        % you may provide other methods required to implement your solver
        % or help the user
        %===============================================================

	   function self = ConstructDomain0(self)

           %update the wave numbers
           id = self.m_DomainArray{1}.m_id;

           start =0;

           if id >0

            intdoms = length(self.m_DomainArray);

           else

               intdoms = length(self.m_DomainArray)-1;

               start =1;

           end

           self.m_DOMS = intdoms;

           totdoms = intdoms+1;

           self.m_kwns_real = zeros(totdoms,1);

           self.m_kwns_imag = zeros(totdoms,1);

           self.m_kwns_real(1) = self.m_kext_real;

           self.m_kwns_imag(1) = self.m_kext_imag;

		   self.m_kwns(1) = self.m_kwns_real(1) + 1i*self.m_kwns_imag(1);

           IntDoms = cell(intdoms,1);

           for i=1:intdoms

               self.m_kwns_real(i+1) = real(self.m_DomainArray{i+start}.GetKwave());

               self.m_kwns_imag(i+1) = imag(self.m_DomainArray{i+start}.GetKwave());

		       self.m_kwns(i+1) = self.m_kwns_real(i+1) + 1i*self.m_kwns_imag(i+1);

               IntDoms{i} = self.m_DomainArray{i+start};

           end

           self.m_DomainArray = cell(totdoms,1);

           for i=1:intdoms

               self.m_DomainArray{i+1} = IntDoms{i};

           end

           InterfzM = cell(1,1);

           count=1;

           for i=1:intdoms
               Nitz = IntDoms{i}.m_NumInfefaces;

               for j =1:Nitz

                   inter = IntDoms{i}.m_Interfaces{j};

                   if((inter.m_exterior==1)&&((inter.m_direction==1)))

                       InterfzM{count} = SetDirection(inter,-1);

                       count = count +1;

                   end

               end

           end

           d0 = Domain(0,self.m_kext_real, self.m_kext_imag,InterfzM);

%            % SQUARE
%            InterfzM{1} = SetDirection(IntDoms{1}.m_Interfaces{1}, -1);
%            InterfzM{2} = SetDirection(IntDoms{4}.m_Interfaces{1}, -1);
%            InterfzM{3} = SetDirection(IntDoms{4}.m_Interfaces{4}, -1);
%            InterfzM{4} = SetDirection(IntDoms{3}.m_Interfaces{4}, -1);
%            InterfzM{5} = SetDirection(IntDoms{3}.m_Interfaces{3}, -1);
%            InterfzM{6} = SetDirection(IntDoms{2}.m_Interfaces{3}, -1);
%            InterfzM{7} = SetDirection(IntDoms{2}.m_Interfaces{2}, -1);
%            InterfzM{8} = SetDirection(IntDoms{1}.m_Interfaces{2}, -1);
%            d0 = Domain(0,self.m_kext_real, self.m_kext_imag,InterfzM);


           self.m_DomainArray{1} = d0;

%            self.m_IsSetup =0;

       	end

       function self = AddDomain(self,dom)

           self.m_DomainArray{end+1} = dom;

           ConstructDomain0(self);
       end


       function self = setKext(self, kreal, kimag)
          self.m_kext_real = kreal;
          self.m_kext_imag = kimag;
          self.m_kext = kreal + 1i*kimag;

       end
       function Display(self)

           hold on;

           for i =2:length(self.m_DomainArray)

               self.m_DomainArray{i}.Display();

           end

           hold off;

        end

    end % end methods

end
