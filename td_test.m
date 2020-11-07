
classdef td_test < incident
    
    properties
        T
        Nt
        m
        tlag
        rk_step
        RK
        dirx
        diry
        wavespeed
        wavespeed2
        setup
        vals
        dxvals
        dyvals
    end
    
    methods

        %-------------------------------------------------
        % constructor
        %-------------------------------------------------
        
        function self = td_test(T, M, m,tlag,dirx, diry,wavespeed, wavespeed2, rk_step, RK)


            % set properties
            self.T = T;
            self.Nt = M;
            self.m = m;
            self.tlag = tlag;
            self.dirx = dirx;
            self.diry = diry;
            self.wavespeed = wavespeed;
            self.wavespeed2 = wavespeed2;
            self.rk_step = rk_step;
            self.RK = RK;
            self.setup = 0;
            
        end
        
        
        function e = POU(self,t,t0,t1)
            % Function used to construct the partition of the 
            % unity. Based on 
            % "Surface scattering in three dimensions: an accelerated
            % high-oder solver" O. P. Bruno and L. A. Kunyansky
            e = zeros(size(t));
            for n=1:numel(t)
                if abs(t(n))<=t0
                    e(n) = 1;
                elseif t0<abs(t(n)) && abs(t(n))<t1
                    x = (abs(t(n))-t0)/(t1-t0); 
                    e(n) = exp(2*exp(-1/x)/(x-1));
                elseif abs(t(n))>=t1
                    e(n) = 0;
                end
            end

        end
        
        function de = dPOU(self,t,t0,t1)
            % Function used to construct the partition of the 
            % unity. Based on 
            % "Surface scattering in three dimensions: an accelerated
            % high-oder solver" O. P. Bruno and L. A. Kunyansky

            de = zeros(size(t));

            for n=1:numel(t)
                if abs(t(n))<=t0
                    %e(n) = 1;
                   de(n) = 0;
                elseif t0<abs(t(n)) && abs(t(n))<t1
                    x = (abs(t(n))-t0)/(t1-t0); 
                    %e(n) = exp(2*exp(-1/x)/(x-1));
                   de(n) = (-2*exp(-1/x+2*exp(-1/x)/(x-1))*(x^2-x+1)/(x^2*(x-1)^2))/(t1-t0)*sign(t(n));
                elseif abs(t(n))>=t1
            %         e(n) = 0;
                   de(n) = 0;
                end
            end

        end    
        
        %-----------------------------------------------
        % evaluate method
        %-----------------------------------------------

        function val = evaluate(self,points)
            
            if self.setup == 0
            
            x = real(points);
            y = imag(points);
            
%             epsil = 10.;
            arg = @(t) self.wavespeed * t - self.tlag  - (x * self.dirx + y * self.diry);

            arg2 = @(t) self.wavespeed2 * t - self.tlag  - (x * self.dirx + y * self.diry);
%             heavyside = @(t) 1./(1+ exp(-epsil * arg(t)));
%             
%             heavyside_bw = @(t) 1./(1+ exp(-epsil * (1-arg(t))));
% 
% 
%             window = @(t) heavyside(t) .* heavyside_bw(t);
% 
% 
%             signal = @(t) sin(16* arg(t)) .* window(t);

            sigma = 2;            
            w = 8;
%             signal = @(t) sin(w*arg(t)) .* exp(-arg(t).^2 / sigma^2) .* (arg(t) > 0);
% 
%             signal = @(t) -0*sin(arg(t)) ./ (1 + exp(-beta * arg(t))) + ...
%                          +sin(arg2(t)) ./ (1 + exp(-beta * arg2(t)));
% 
% 
%             signal = @(t) -0*sin(arg(t)) .* (1-self.POU(arg(t),0.2, 2)) .* (arg(t)>0) + ...
%                          +sin(arg2(t)) .* (1-self.POU(arg2(t), 0.2, 2)) .* (arg2(t)>0);
                     
                     


            
%             signal = @(t) sin(w*arg(t)) .*(1-self.POU(arg(t)-3,0.2, 2)) .* (arg(t) > 0);

            signal = @(t) sin(w*arg(t)) .*(1- self.POU(arg(t),0.2, 5)) .* (arg(t) > 0);
            % intialize return array
            size_val = size(points);
            size_val(length(size_val) + 1) = self.Nt+1;
            
            vals=zeros(size_val);
  

            % compute field
%             time = linspace(0, self.T, self.Nt);
            dt = self.T/self.Nt;
            time = 0:dt:self.T;
            R = eps^(0.5/(self.Nt+1));
            for n = 1:self.Nt+1
                
                vals(:, :, n) = R^(n-1) * signal(time(n) + self.RK*self.rk_step);
            
            end

            
             
            vals = fft(vals, [], length(size_val));
            
            val = vals(:, :, self.m);
            
            self.vals = vals;
            
%             self.setup = 1;
             
           
            
            
            else
                
                
           
            
            val = self.vals(:, :, self.m);
            
            end
            
        end
        
        %-----------------------------------------------
        % evaluate method
        %-----------------------------------------------

        function [dxm,dym] = evaluateGradient(self,points)
            
            
            if self.setup == 0
            
            x = real(points);
            y = imag(points);
            
%             epsil = 10.;
            arg = @(t) self.wavespeed * (t - self.tlag)  - (x * self.dirx + y * self.diry);
            
            arg2 = @(t) self.wavespeed2 * t - self.tlag  - (x * self.dirx + y * self.diry);
% 
%             heavyside = @(t) 1./(1+ exp(-epsil * arg(t)));
%             
%             heavyside_bw = @(t) 1./(1+ exp(-epsil * (1-arg(t))));
% 
%             dheavyside = @(t) 1./(1+exp(-arg(t)*epsil)).^2*epsil.*exp(-arg(t)*epsil);
% 
%             dheavyside_bw = @(t) 1./(1+exp(-(1-arg(t))*epsil)).^2*epsil.*exp(-(1-arg(t))*epsil);
% 
%             window = @(t) heavyside(t) .* heavyside_bw(t);
% 
%             windowp = @(t) dheavyside(t) .* heavyside_bw(t) - ...
%                            heavyside(t) .* dheavyside_bw(t);
% 

%             signalp = @(t) 16*cos(16* arg(t)) .*window(t) + sin(16* arg(t)) .* windowp(t);
%             w = 1;
%             sigma = 2;
% %             signalp = @(t) w*cos(w* arg(t)) .* exp(-arg(t).^2 / sigma^2) .* (arg(t) > 0)- 2 * (arg(t)) .* sin(w* arg(t)) .* exp(- (arg(t)).^2 / sigma^2) / sigma^2.* (arg(t) > 0) ;
%             beta = 10;
% %             
% %             signalp = @(t) 0*cos(arg(t)) ./ (1+ exp(-beta*arg(t))) + ...
% %                0*sin(arg(t))./(1+exp(-arg(t)*beta)).^2*beta.*exp(-arg(t)*beta) +...
% %                cos(arg2(t)) ./ (1+ exp(-beta*arg2(t))) + ...
% %                sin(arg2(t))./(1+exp(-arg2(t)*beta)).^2*beta.*exp(-arg2(t)*beta);
% 

%             w = 4;
%             signalp = @(t) -0*cos(arg(t)) .* (1-self.POU(arg(t), 0.2, 2)) .*(arg(t)>0) + ...
%                -0*sin(arg(t)).* (-self.dPOU(arg(t), 0.2, 2)) .*(arg(t)>0) +...
%                cos(arg2(t)) .* (1-self.POU(arg2(t), 0.2, 2)) .*(arg2(t)>0) + ...
%                sin(arg2(t)) .* (-self.dPOU(arg2(t), 0.2, 2)) .*(arg2(t)>0);
% %            
%            

            w = 8;
            signalp = @(t) w*cos(w*arg(t)) .* (1-self.POU(arg(t), 0.2, 5)) .* (arg(t)> 0) + ...
               -sin(w*arg(t)) .* self.dPOU(arg(t), 0.2, 5) .* (arg(t)>0);

%             signalp = @(t) w*cos(w*arg(t)) .* self.POU(arg(t)-4, 0.2, 4) .* (arg(t)> 0) + ...
%                -sin(w*arg(t)) .* self.dPOU(arg(t)-4, 0.2, 4) .* (arg(t)>0);
% intialize return array
            size_val = size(points);
            size_val(length(size_val) + 1) = self.Nt+1;
            
            % intialize return array
            dx=zeros(size_val);
            dy=zeros(size_val);



            
%             time = linspace(0, self.T, self.Nt);
            dt = self.T/self.Nt;
            time = 0:dt:self.T;
            R = eps^(0.5/(self.Nt+1));
            for n = 1:self.Nt+1
                
                % now compute gradient


                dx(:, :, n) = - self.dirx * R^(n-1) * signalp(time(n)+ self.RK*self.rk_step);
                dy(:, :, n) = - self.diry * R^(n-1) * signalp(time(n)+ self.RK*self.rk_step);

            end
            

            
            dx = fft(dx, [], length(size_val));
            
            
            dy = fft(dy, [], length(size_val));
            
            self.dxvals = dx;
            self.dyvals = dy;
            
            dxm = dx(:, :, self.m);
            
            dym = dy(:, :, self.m);
            
            self.setup = 1;
            
            
            
            else
                
    
            
            dxm = self.dxvals(:, :, self.m);
            
            dym = self.dyvals(:, :, self.m);
                
                
                
            end
            
        end
        
        
        function self = Setm(self, m)
           
            self.m = m;
            
        end
        
        
    end
        % end methods
    
end
