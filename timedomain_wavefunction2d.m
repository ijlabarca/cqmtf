classdef timedomain_wavefunction2d < incident

    properties
        time
        tlag
        dirx
        diry
        wavespeed
        m_time
    end

    methods

        %-------------------------------------------------
        % constructor
        %-------------------------------------------------

        function self = timedomain_wavefunction2d(time,tlag,dirx, diry,wavespeed)


            % set properties
            self.time = time;
            self.tlag = tlag;
            self.dirx = dirx;
            self.diry = diry;
            self.wavespeed = wavespeed;
            self.m_time = 1;

        end


        %-----------------------------------------------
        % smooth version of heavyside function
        %-----------------------------------------------
%
%         function val = heavyside(self, points, rev)
%
%
%             x = real(points);
%             y = imag(points);
%             epsil = 10;
%
%             arg = self.wavespeed * (self.time - self.tlag)  - (x * self.dirx + y * self.diry);
%
%             if rev == 1
%
%                arg = 1- arg;
%             end
%
%             val = 1./(1+exp(-arg*epsil));
%
%         end

        %-----------------------------------------------
        % derivative of heavyside function
        %-----------------------------------------------

%         function val = dheavyside(self,points, rev)
%
%             x = real(points);
%             y = imag(points);
%             epsil = 10;
%             arg = self.wavespeed * (self.time - self.tlag)  - (x * self.dirx + y * self.diry);
%
%             if rev == 1
%
%                arg = 1- arg;
%
%             end
%
%             val = 1./(1+exp(-arg*epsil)).^2*epsil.*exp(-arg*epsil);
%
%         end



        %-----------------------------------------------
        % window function
        %-----------------------------------------------

%         function val = window(self, points)
%
%             val = self.heavyside(points, 0) .* self.heavyside(points, 1);
%
% %             val = self.heavyside(points, 0);
%
%         end


        %-----------------------------------------------
        % derivative of window function
        %-----------------------------------------------

%         function val = windowp(self, points)
%
%             val = self.dheavyside(points, 0).*self.heavyside(points, 1) ...
%                  -self.heavyside(points, 0).*self.dheavyside(points, 1);
%
%
% %             val = self.dheavyside(points, 0);
%
%         end


        %-----------------------------------------------
        % sine
        %-----------------------------------------------

%         function val = sine(self, points)
%
%
%             x = real(points);
%             y = imag(points);
%
%             arg = self.wavespeed * (self.time - self.tlag)  - (x * self.dirx + y * self.diry);
%
%
%             val = sin(16* arg);
%
%         end

        %-----------------------------------------------
        % derivative of sine
        %-----------------------------------------------

%         function val = dsine(self,points)
%
%             x = real(points);
%             y = imag(points);
%             arg = self.wavespeed * (self.time - self.tlag)  - (x * self.dirx + y * self.diry);
%
%             val = 16*cos(16*arg);
%
%         end
%

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

        function val = signal(self, points)
            x = real(points);
            y = imag(points);
            arg = self.wavespeed * (self.time - self.tlag)  - (x * self.dirx + y * self.diry);
            sigma = 2;
            w = 8;
            val = sin(w*arg) .* (1-self.POU(arg, 0.2, 5)) .* (arg > 0);

%             val = sin(w*arg) .* self.POU(arg-4, 0.2, 4) .* (arg > 0);
%             val = sin(w*arg) .* exp(-arg.^2 / sigma^2) .* (arg > 0);

        end


        function val = signalp(self, points)

            x = real(points);
            y = imag(points);
            arg = self.wavespeed * (self.time - self.tlag)  - (x * self.dirx + y * self.diry);
            w = 8;
            sigma = 2*sqrt(2/3);
            val = w*cos(w* arg) .* exp(-arg.^2 / sigma^2) - 2 * (arg) .* sin(w* arg) .* exp(- (arg).^2 / sigma^2) / sigma^2 ;


        end


        %-----------------------------------------------
        % evaluate method
        %-----------------------------------------------

        function val = evaluate(self,points,mask)

            % intialize return array
            val=zeros(size(points));
%             val = NaN(size(points));

            % apply mask if necessary
            if nargin>2
                points=points(mask);
            end


            % compute field
%             v = self.sine(points).* self.window(points);

            v = self.signal(points);

            % insert values into the return array
            if nargin>2
                val(mask)=v;
            else
                val=v;
            end

        end

        %-----------------------------------------------
        % evaluate method
        %-----------------------------------------------

        function [dx,dy] = evaluateGradient(self,points,mask)

            % intialize return array
            dx=zeros(size(points));
            dy=zeros(size(points));

            % apply mask if necessary
            if nargin>2
                points=points(mask);
            end

%             df = self.sine(points).* self.windowp(points)+ self.dsine(points).* self.window(points);

            df = self.signalp(points);

            % now compute gradient
            vx = - self.dirx * df;
            vy = - self.diry * df;

            % insert values into the return array
            if nargin>2
                dx(mask) = vx;
                dy(mask) = vy;
            else
                dx = vx;
                dy = vy;
            end

        end

    end % end methods

end
