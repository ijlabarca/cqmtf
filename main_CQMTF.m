clear all
close all
clc

%%
%-----------------------------------------
% set main parameters
%-----------------------------------------
eme = 0;
ene = 0;


M = 200; % timesteps
disp('-----------------------------------------------')
% wavespeeds
Nmax = 20;
% M = emes(eme);

c_0 = 1.;
c_1 = 0.5;
c_2 = 0.25;

% origin
center = 0;

% domain wiht one elipse.
a=1;
b=1;



%left half of the elipse;
x =@(t)a*cos(pi/2*(t+1)+pi-pi/2);
y = @(t) b*sin(pi/2*(t+1)+pi-pi/2);
dx = @(t) -a*pi/2*sin(pi/2*(t+1)+pi-pi/2);
dy =@(t) b*pi/2*cos(pi/2*(t+1)+pi-pi/2);
ddx = @(t) -a*pi*pi/4*cos(pi/2*(t+1)+pi-pi/2);
ddy = @(t) -b*pi*pi/4*sin(pi/2*(t+1)+pi-pi/2);

%construction of the curve
c2 = Interface(x,y,dx,dy,ddx,ddy);

%right half of the elipse;
x =@(t)a*cos(pi/2*(t+1)-pi/2);
y = @(t) b*sin(pi/2*(t+1)-pi/2);
dx = @(t) -a*pi/2*sin(pi/2*(t+1)-pi/2);
dy =@(t) b*pi/2*cos(pi/2*(t+1)-pi/2);
ddx = @(t) -a*pi*pi/4*cos(pi/2*(t+1)-pi/2);
ddy = @(t) -b*pi*pi/4*sin(pi/2*(t+1)-pi/2);
c1 = Interface(x,y,dx,dy,ddx,ddy);

%separation between the two halves:
x=@(t) 0;
y= @(t) -b*t;
dx =@(t) 0;
dy =@(t) -b*1;
ddx =@(t) 0;
ddy =@(t) 0;
% we add a 0 at the end to indicate that this curve is not on the exterior
% domain
c3 =Interface(x,y,dx,dy,ddx,ddy,0);

%separtion oriented in oposite direction.
c4 = c3.SetDirection(-1);


%first domain Left domain.
%NOT USE id =0!!!!
 Interfz = cell(2,1);
 Interfz{1} =c2;
 Interfz{2} =c4;
 k1=4;
 d1 =  Domain(1,0,-1.0*k1,Interfz);



 %Second domain Right domain.
 Interfz = cell(2,1);
 Interfz{1} =c3;
 Interfz{2} =c1;
 k2 =8;
 d2 =  Domain(2,0,-1.0*k2,Interfz);

 % We save all the domains in one cell array
 Domains = cell(2,1);
 Domains{1} =d1;
 Domains{2} =d2;

%     M = 100;
T = 10;
dt = T/M;
time = 0:dt:T;
% time = linspace(0, T, M);
tlag = 3.5;
dirx = 0;
diry = -1;
% BDF2 
p = @(z) 1.5-2*z+0.5*z.^2;

omega = exp(2*pi*1i/(M+1));
R = eps^(0.5/(M+1));

% parfor m = 0:M
% %     inc{m + 1} = timedomain_wavefunction2d(time(m +1),tlag,dirx, diry, c_0);
% 
%     inc{m + 1} = tdwavefunction2d(T, M, m+1,tlag,dirx, diry, c_0);
% end



inc{1} = tdwavefunction2d(T, M, 1,tlag,dirx, diry, c_0);

% setup solver with empty incident field for now
kwave_real = 0;
kwave_imag = -4*pi;

wavenums = p(R*omega.^(-(0:floor(M/2))))/ (c_0 * T / M);
maxK = max(abs(wavenums));
% Nmax = ceil(1.1*maxK)+7
%%

solver = solverMtf28(c_0,kwave_real,kwave_imag,inc,Domains);


        tic
solver.setupb(Nmax, M)

disp('Lado derecho')
        toc
% solver.m_b = bsxfun(@times,solver.m_b,R.^(0:M-1));  % scaling
% solver.m_b = fft(solver.m_b, [], 2);    % dft by columns
d = size(solver.m_b);
g = zeros(d(1),M+1);

%%
disp(' ')
disp('TDMTF-solve');
tic
hbar = parfor_progressbar(M/2+1,'Solving Linear Systems...'); %create the progress bar 
for l=0:M/2       %compute half the sequence

%             disp('-------------------------------------------------');

%     if l == 0
% 
%         tic
% 
% 
%     end
%     tic
        


    if norm(solver.m_b(:,l+1)) > 1e-12
    k0 = p(R*omega^(-l))/(c_0 * T/M);
    k1 = p(R*omega^(-l))/(c_1 * T/M);
    k2 = p(R*omega^(-l))/(c_2 * T/M);

    %Set wave numbers for each subdomain
    solver.m_kext_real = real(k0);

    solver.m_kext_imag = imag(k0);

    solver.m_kext = k0;


    solver.m_DomainArray{1}.m_KwaveNumber = k0;
    solver.m_DomainArray{2}.m_KwaveNumber = k1;
    solver.m_DomainArray{3}.m_KwaveNumber = k2;
    solver.ConstructDomain0();

    solver.setup(Nmax)

    g(:,l+1)= solver.m_A \ solver.m_b(:,l+1);  
    

%     disp(solver.CalderonError(g(:, l+1)) );


%     if l == 0
%        toc 
%     end

%     [g(:,l+1), ~] = gmres(solver.m_A, solver.m_b(:,l+1));
% 
%             disp(k0);
%     solver.JumpError(g(:, l+1), l, c_0, T, M, tlag, dirx, diry, d(1))

%     norm(g(:,l+1))

    end
    
    hbar.iterate(1);
toc
end
toc
close(hbar);
disp(['Number of timesteps: ', num2str(M)]);

disp(' ')

g(:,M+2-(1:floor(M/2)))=conj(g(:,2:floor(M/2)+1));  % mirror the hermitian seq
g=real(ifft(g,[],2));                            % w
g=bsxfun(@times,g,R.^(-(0:M)));                % g




%%

disp('Computing Traces');
% Ng = 2*max(solver.m_N{1});

Ng = 101;
[xg,wg] = lgwt(Ng,-1,1);
Un =@(n,t) sin((n+1)*acos(t))./sin(acos(t));

tD = zeros(3,2, Ng, M+1);
tN = zeros(3,2, Ng, M+1);

deltat = T/M;

uinc = zeros(2,Ng, M+1);
dnuinc = zeros(2,Ng, M+1);


lambdaD = zeros(3, size(g, 1)/6, M+1);
lambdaN = zeros(3, size(g, 1)/6, M+1);


lambdaD(1, :, :) = g(1:size(g, 1)/6, :);
lambdaN(1, :, :) = g(size(g, 1)/6 + 1:size(g, 1)/3, :);

lambdaD(2, :, :) = g(size(g, 1)/3 + 1:size(g, 1)/2, :);
lambdaN(2, :, :) = g(size(g, 1)/2 + 1:2*size(g, 1)/3, :);

lambdaD(3, :, :) = g(2*size(g, 1)/3 + 1:5*size(g, 1)/6, :);
lambdaN(3, :, :) = g(5*size(g, 1)/6 + 1:size(g, 1), :);
%%
tic
for dom = 1:3


    counter = 0;

     NIntrfz0 = solver.m_DomainArray{dom}.m_NumInfefaces;

     for jj=1:NIntrfz0



               [nx,ny] = solver.m_DomainArray{dom}.normal(xg,jj);

               [x,y] = solver.m_DomainArray{dom}.geo(xg,jj);

               z = x + 1i*y;

               J = solver.m_DomainArray{dom}.J(xg, jj);



               for m = 1:M+1


                   for n = 0:2*solver.m_N{dom}(jj)


                        for gg=1:Ng



                           if length(J) == 1


                               tD(dom,jj, gg, m) = tD(dom,jj, gg, m) + (lambdaD(dom, counter + n + 1, m) * Un(n, xg(gg))./J);  

                               tN(dom,jj, gg, m) = tN(dom,jj, gg, m) + (lambdaN(dom, counter + n + 1, m) * Un(n, xg(gg))./J); 

                           else 

                               tD(dom,jj, gg, m) = tD(dom,jj, gg, m) + (lambdaD(dom, counter + n + 1, m) * Un(n, xg(gg))./J(gg));  

                               tN(dom,jj, gg, m) = tN(dom,jj, gg, m) + (lambdaN(dom, counter + n + 1, m) * Un(n, xg(gg))./J(gg)); 


                           end
                        end


                   end


                   if dom == 1
                       onda = timedomain_wavefunction2d((m-1)*deltat, tlag, dirx, diry, c_0);
                       [dx, dy] = onda.evaluateGradient(z);

                       dnuinc(jj,:,m) = nx .* dx + ny .* dy;
                       uinc(jj,:, m) = onda.evaluate(z);

                   end


               end

               counter = counter + 2*solver.m_N{dom}(jj) + 1;

     end

end

toc

plot_fields;


%% % Plot Traces
% % 
% plotTraces(deltat, tD,0,uinc, 1, 2, 1, 1);
% 
% 
% plotTraces(deltat, tD,0,uinc, 2, 3, 2, 1);
% 
% 
% plotTraces(deltat, tD,0,uinc, 1, 3, 2, 2);


