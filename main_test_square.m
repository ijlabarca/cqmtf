
clear 
close all
clc

%%
%-----------------------------------------
% set main parameters
%-----------------------------------------
eme = 0;
ene = 0;


% M = 10; % timesteps


emes =[50 100 200 400];
Mmax = 800;
% for Nmax = [20, eNes]


for M = [Mmax emes]


% M = 1280;

disp('-----------------------------------------------')
% wavespeeds
Nmax = 8;
% M = emes(eme);



c_0 = 1.;
c_1 = 0.5;
c_2 = 0.25;
c_3 = 0.5;
c_4 = 0.25;
% origin
center = 0;

% domain wiht one elipse.
a=0.5;
b=0.5;



%% c1:
x=@(t) a/2*(1-t);
y= @(t) b;
dx =@(t) -a/2;
dy =@(t) 0;
ddx =@(t) 0;
ddy =@(t) 0;
c1 = Interface(x,y,dx,dy,ddx,ddy);

%% c2:
x=@(t) -a/2*(1+t);
y= @(t) b;
dx =@(t) -a/2;
dy =@(t) 0;
ddx =@(t) 0;
ddy =@(t) 0;
c2 = Interface(x,y,dx,dy,ddx,ddy);

%% c3:
x=@(t) -a;
y= @(t) b/2*(1-t);
dx =@(t) 0;
dy =@(t) -b/2;
ddx =@(t) 0;
ddy =@(t) 0;
c3 = Interface(x,y,dx,dy,ddx,ddy);


%% c4:
x=@(t) -a;
y= @(t) -b/2*(1+t);
dx =@(t) 0;
dy =@(t) -b/2;
ddx =@(t) 0;
ddy =@(t) 0;
c4 = Interface(x,y,dx,dy,ddx,ddy);


%% c5:
x=@(t) -a/2*(1-t);
y= @(t) -b;
dx =@(t) a/2;
dy =@(t) 0;
ddx =@(t) 0;
ddy =@(t) 0;
c5 = Interface(x,y,dx,dy,ddx,ddy);

%% c6:
x=@(t) a/2*(1+t);
y= @(t) -b;
dx =@(t) a/2;
dy =@(t) 0;
ddx =@(t) 0;
ddy =@(t) 0;
c6 = Interface(x,y,dx,dy,ddx,ddy);

%% c7:
x=@(t) a;
y= @(t) -b/2*(1-t);
dx =@(t) 0;
dy =@(t) b/2;
ddx =@(t) 0;
ddy =@(t) 0;
c7 = Interface(x,y,dx,dy,ddx,ddy);

%% c8:
x=@(t) a;
y= @(t) b/2*(1+t);
dx =@(t) 0;
dy =@(t) b/2;
ddx =@(t) 0;
ddy =@(t) 0;
c8 = Interface(x,y,dx,dy,ddx,ddy);

%% c9:
x=@(t) a/2*(1-t);
y= @(t) 0;
dx =@(t) -a/2;
dy =@(t) 0;
ddx =@(t) 0;
ddy =@(t) 0;
c9 = Interface(x,y,dx,dy,ddx,ddy, 0);

%% c10:
x=@(t) -a/2*(1+t);
y= @(t) 0;
dx =@(t) -a/2;
dy =@(t) 0;
ddx =@(t) 0;
ddy =@(t) 0;
c10 = Interface(x,y,dx,dy,ddx,ddy, 0);



%% c11:
x=@(t) 0;
y= @(t) b/2*(1-t);
dx =@(t) 0;
dy =@(t) -b/2;
ddx =@(t) 0;
ddy =@(t) 0;
c11 = Interface(x,y,dx,dy,ddx,ddy, 0);


%% c12:
x=@(t) 0;
y= @(t) -b/2*(1+t);
dx =@(t) 0;
dy =@(t) -b/2;
ddx =@(t) 0;
ddy =@(t) 0;
c12 = Interface(x,y,dx,dy,ddx,ddy, 0);


% top left.
%NOT USE id =0!!!!
 Interfz = cell(4,1);
 Interfz{1} =c2;
 Interfz{2} =c3;
 Interfz{3} =c10.SetDirection(-1);
 Interfz{4} =c11.SetDirection(-1);
 k1=4;
 d1 =  Domain(1,0,-1.0*k1,Interfz);


% bottom left.
%NOT USE id =0!!!!
 Interfz = cell(4,1);
 Interfz{1} =c10;
 Interfz{2} =c4;
 Interfz{3} =c5;
 Interfz{4} =c12.SetDirection(-1);
 k1=4;
 d2 =  Domain(2,0,-1.0*k1,Interfz);

% bottom right.
%NOT USE id =0!!!!
 Interfz = cell(4,1);
 Interfz{1} =c9;
 Interfz{2} =c12;
 Interfz{3} =c6;
 Interfz{4} =c7;
 k1=4;
 d3 =  Domain(3,0,-1.0*k1,Interfz);

% top right.
%NOT USE id =0!!!!
 Interfz = cell(4,1);
 Interfz{1} =c1;
 Interfz{2} =c11;
 Interfz{3} =c9.SetDirection(-1);
 Interfz{4} =c8;
 k1=4;
 d4 =  Domain(4,0,-1.0*k1,Interfz);


 % We save all the domains in one cell array
 Domains = cell(4,1);
 Domains{1} =d1;
 Domains{2} =d2;
 Domains{3} =d3;
 Domains{4} =d4;

%     M = 100;
T = 10; % 10
dt = T/M;
time = 0:dt:T;
% time = linspace(0, T, M);
tlag = 0.5;
dirx = sqrt(0.5);
diry = -sqrt(0.5);
% BDF2 
p = @(z) 1.5-2*z+0.5*z.^2;

omega = exp(2*pi*1i/(M+1));
R = eps^(0.5/(M+1));

inc{1} = td_test(T, M, 1,tlag,dirx, diry, c_0, c_1, 0, 0);
%inc2{1} = td_test(T, M, 1,tlag,sqrt(0.5), sqrt(0.5), c_1);

% for m = 0:M
% %     inc{m + 1} = timedomain_wavefunction2d(time(m +1),tlag,dirx, diry, c_0);
% 
%     inc{m + 1} = tdwavefunction2d(T, M, m+1,tlag,dirx, diry, c_0);

% end
% setup solver with empty incident field for now


kwave_real = 0;
kwave_imag = -4*pi;

wavenums = p(R*omega.^(-(0:floor(M/2))))/ (c_0 * T / M);
maxK = max(abs(wavenums));
% Nmax = ceil(1.1*maxK)+7
%%

solver = solverMtf28(c_0,kwave_real,kwave_imag,inc,Domains);

%solver2 = solverMtf28(c_0,kwave_real,kwave_imag,inc2,Domains);

        tic
solver.setupb(Nmax, M)
%solver2.setupb(Nmax, M)
disp('Lado derecho')
        toc
% solver.m_b = bsxfun(@times,solver.m_b,R.^(0:M-1));  % scaling
% solver.m_b = fft(solver.m_b, [], 2);    % dft by columns
d = size(solver.m_b);
g = zeros(d(1),M+1);

% 
% step = size(solver2.m_b,1)/6;
% start = size(solver2.m_b,1)/6;
% 
% solver2.m_b(start+1:start+step, :) = -solver2.m_b(start+1:start+step, :); 
% start = start + 2*step;
% 
% solver2.m_b(start+1:start+step, :) = -solver2.m_b(start+1:start+step, :); 
% start = start + 2*step;
% 
% solver2.m_b(start+1:start+step, :) = -solver2.m_b(start+1:start+step, :); 
% cald_error = zeros(size(emes, 1)+1, max(emes));

%%
disp(' ')
disp('TDMTF-solve');
tic
% hbar = parfor_progressbar(M/2+1,'Solving Linear Systems...'); %create the progress bar 
for l=0:M/2       %compute half the sequence

%             disp('-------------------------------------------------');

%     if l == 0
% 
%         tic
% 
% 
%     end
%     tic
        


    if norm(solver.m_b(:,l+1)) > 1e-14
    k0 = p(R*omega^(-l))/(c_0 * T/M);
    k1 = p(R*omega^(-l))/(c_1 * T/M);
    k2 = p(R*omega^(-l))/(c_2 * T/M);
    k3 = p(R*omega^(-l))/(c_3 * T/M);
    k4 = p(R*omega^(-l))/(c_4 * T/M);
    %Set wave numbers for each subdomain
    solver.m_kext_real = real(k0);

    solver.m_kext_imag = imag(k0);

    solver.m_kext = k0;


    solver.m_DomainArray{1}.m_KwaveNumber = k0;
    solver.m_DomainArray{2}.m_KwaveNumber = k1;
    solver.m_DomainArray{3}.m_KwaveNumber = k2;
    solver.m_DomainArray{4}.m_KwaveNumber = k3;
    solver.m_DomainArray{5}.m_KwaveNumber = k4;
    solver.ConstructDomain0();

    solver.setup(Nmax)
%     rhs = solver.m_b(:,l+1) - solver2.m_b(:,l+1);
    rhs = solver.m_b(:,l+1);
    
    g(:,l+1)= solver.m_A \ rhs;  
    
%     if eme ~= 0
%     cald_error(eme+1, l+1) = solver.CalderonError(g(:, l+1));
    disp(l);
%     clear mex_funs
%     end

%     if l == 0
%        toc 
%     end

%     [g(:,l+1), ~] = gmres(solver.m_A, solver.m_b(:,l+1));
% 
%             disp(k0);
%     solver.JumpError(g(:, l+1), l, c_0, T, M, tlag, dirx, diry, d(1))

%     norm(g(:,l+1))

    end
    
%     hbar.iterate(1);
toc
end
% toc
% close(hbar);
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
% xg = linspace(-1, 1, Ng+2);
% xg = xg(2:end-1).';
Un =@(n,t) sin((n+1)*acos(t))./sin(acos(t));

tD0 = zeros(1,8, Ng, M+1);
tN0 = zeros(1,8, Ng, M+1);


tD = zeros(4,4, Ng, M+1);
tN = zeros(4,4, Ng, M+1);

deltat = T/M;
dim = size(g, 1)/12;

uinc = zeros(8,Ng, M+1);
dnuinc = zeros(8,Ng, M+1);


lambdaD0 = zeros(1, 2*dim, M+1);
lambdaN0 = zeros(1, 2*dim, M+1);


lambdaD = zeros(4, dim, M+1);
lambdaN = zeros(4, dim, M+1);

gg = g;

lambdaD0(1, :, :) = gg(1:2*dim, :); gg(1:2*dim, :) = [];
lambdaN0(1, :, :) = gg(1:2*dim, :); gg(1:2*dim, :) = [];


lambdaD(1, :, :) = gg(1:dim, :); gg(1:dim, :) = [];
lambdaN(1, :, :) = gg(1:dim, :); gg(1:dim, :) = [];


lambdaD(2, :, :) = gg(1:dim, :); gg(1:dim, :) = [];
lambdaN(2, :, :) = gg(1:dim, :); gg(1:dim, :) = [];

lambdaD(3, :, :) = gg(1:dim, :); gg(1:dim, :) = [];
lambdaN(3, :, :) = gg(1:dim, :); gg(1:dim, :) = [];

lambdaD(4, :, :) = gg(1:dim, :); gg(1:dim, :) = [];
lambdaN(4, :, :) = gg(1:dim, :); gg(1:dim, :) = []

clear gg
%%
tic
for dom = 1:5


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


                        if dom == 1
                            
                           if length(J) == 1


                               tD0(dom,jj, gg, m) = tD0(dom,jj, gg, m) + (lambdaD0(dom, counter + n + 1, m) * Un(n, xg(gg))./J);  

                               tN0(dom,jj, gg, m) = tN0(dom,jj, gg, m) + (lambdaN0(dom, counter + n + 1, m) * Un(n, xg(gg))./J); 

                           else 

                               tD0(dom,jj, gg, m) = tD0(dom,jj, gg, m) + (lambdaD0(dom, counter + n + 1, m) * Un(n, xg(gg))./J(gg));  

                               tN0(dom,jj, gg, m) = tN0(dom,jj, gg, m) + (lambdaN0(dom, counter + n + 1, m) * Un(n, xg(gg))./J(gg)); 


                           end
                           
                        else
                            
                           if length(J) == 1


                               tD(dom-1,jj, gg, m) = tD(dom-1,jj, gg, m) + (lambdaD(dom-1, counter + n + 1, m) * Un(n, xg(gg))./J);  

                               tN(dom-1,jj, gg, m) = tN(dom-1,jj, gg, m) + (lambdaN(dom-1, counter + n + 1, m) * Un(n, xg(gg))./J); 

                           else 

                               tD(dom-1,jj, gg, m) = tD(dom-1,jj, gg, m) + (lambdaD(dom-1, counter + n + 1, m) * Un(n, xg(gg))./J(gg));  

                               tN(dom-1,jj, gg, m) = tN(dom-1,jj, gg, m) + (lambdaN(dom-1, counter + n + 1, m) * Un(n, xg(gg))./J(gg)); 


                           end
                           
                        end

                        end


                   end


%                    if dom == 1
%                        onda = timedomain_test((m-1)*deltat, tlag, dirx, diry, c_0);
%                        [dx, dy] = onda.evaluateGradient(z);
% 
%                        tN_ref(1,jj,:,m) = nx .* dx + ny .* dy;
%                        tD_ref(1,jj,:, m) = onda.evaluate(z);
%                        
%                    elseif dom == 2
% 
%                        onda = timedomain_test((m-1)*deltat, tlag, dirx, diry, c_1);
%                        [dx, dy] = onda.evaluateGradient(z);
% 
%                        tN_ref(2,jj,:,m) = nx .* dx + ny .* dy;
%                        tD_ref(2,jj,:, m) = onda.evaluate(z);
%                    elseif dom == 3
%                        
%                        onda = timedomain_test((m-1)*deltat, tlag, dirx, diry, c_2);
%                        [dx, dy] = onda.evaluateGradient(z);
% 
%                        tN_ref(3,jj,:,m) = nx .* dx + ny .* dy;
%                        tD_ref(3,jj,:, m) = onda.evaluate(z);
%                        
%                        
%                    end


               end

               counter = counter + 2*solver.m_N{dom}(jj) + 1;

     end

end

toc
%% Errors

% solver.m_kext_real = 10;
% solver.m_kext_imag = 0;
% solver.m_kext = 10;
% solver.m_DomainArray{1}.m_KwaveNumber = 10;
% solver.m_DomainArray{2}.m_KwaveNumber = 10;
% solver.m_DomainArray{3}.m_KwaveNumber = 10;
% solver.ConstructDomain0();
% solver.setup((Ng-1)/2);
% 
% Ai = solver.m_Ai;

load('Ai.mat')

% %%
% disp('------------------------------------------------');
% disp(Nmax);
% disp(M);
% disp('Error Traza Dirichlet')
% 
% if eme ~= 0
% jumpError_01 = jumpError(Ai, deltat, tD,0, uinc, 1, 2, 1, 1);
% jumpError_02 = jumpError(Ai, deltat, tD,0, uinc, 1, 3, 2, 2);
% jumpError_12 = jumpError(Ai, deltat, tD,0, uinc, 2, 3, 2, 1);
% 
% jumpErrorD(eme, 1) = jumpError_01;
% 
% jumpErrorD(eme, 2) = jumpError_02;
% 
% jumpErrorD(eme, 3) = jumpError_12;
% % 
% % disp('------------------------------------------------');
% % 
% % disp('Error Traza Neumann')
% jumpError_01 = jumpError(Ai, deltat, tN,1, dnuinc, 1, 2, 1, 1);
% jumpError_02 = jumpError(Ai, deltat, tN,1, dnuinc, 1, 3, 2, 2);
% jumpError_12 = jumpError(Ai, deltat, tN,1, dnuinc, 2, 3, 2, 1);
% 
% 
% jumpErrorN(eme, 1) = jumpError_01;
% 
% jumpErrorN(eme, 2) = jumpError_02;
% 
% jumpErrorN(eme, 3) = jumpError_12;
% % 
% end


%plotTraces(deltat, tD,0,tD_ref, 1, 2, 1, 1);
% 
% 
%plotTraces(deltat, tD,0,tD_ref, 2, 3, 2, 1);
% 
% 
%plotTraces(deltat, tD,0,tD_ref, 1, 3, 2, 2);

%return



%% c1:
x=@(t) a/2*(1-t);
y= @(t) b*ones(numel(t), 1);
dx =@(t) -a/2*ones(numel(t), 1);
dy =@(t) 0*ones(numel(t), 1);
ddx =@(t) 0*ones(numel(t), 1);
ddy =@(t) 0*ones(numel(t), 1);
c1 = Interface(x,y,dx,dy,ddx,ddy);

%% c2:
x=@(t) -a/2*(1+t);
y= @(t) b*ones(numel(t), 1);
dx =@(t) -a/2*ones(numel(t), 1);
dy =@(t) 0*ones(numel(t), 1);
ddx =@(t) 0*ones(numel(t), 1);
ddy =@(t) 0*ones(numel(t), 1);
c2 = Interface(x,y,dx,dy,ddx,ddy);

%% c3:
x=@(t) -a*ones(numel(t), 1);
y= @(t) b/2*(1-t);
dx =@(t) 0*ones(numel(t), 1);
dy =@(t) -b/2*ones(numel(t), 1);
ddx =@(t) 0*ones(numel(t), 1);
ddy =@(t) 0*ones(numel(t), 1);
c3 = Interface(x,y,dx,dy,ddx,ddy);


%% c4:
x=@(t) -a*ones(numel(t), 1);
y= @(t) -b/2*(1+t);
dx =@(t) 0*ones(numel(t), 1);
dy =@(t) -b/2*ones(numel(t), 1);
ddx =@(t) 0*ones(numel(t), 1);
ddy =@(t) 0*ones(numel(t), 1);
c4 = Interface(x,y,dx,dy,ddx,ddy);


%% c5:
x=@(t) -a/2*(1-t);
y= @(t) -b*ones(numel(t), 1);
dx =@(t) a/2*ones(numel(t), 1);
dy =@(t) 0*ones(numel(t), 1);
ddx =@(t) 0*ones(numel(t), 1);
ddy =@(t) 0*ones(numel(t), 1);
c5 = Interface(x,y,dx,dy,ddx,ddy);

%% c6:
x=@(t) a/2*(1+t);
y= @(t) -b*ones(numel(t), 1);
dx =@(t) a/2*ones(numel(t), 1);
dy =@(t) 0*ones(numel(t), 1);
ddx =@(t) 0*ones(numel(t), 1);
ddy =@(t) 0*ones(numel(t), 1);
c6 = Interface(x,y,dx,dy,ddx,ddy);

%% c7:
x=@(t) a*ones(numel(t), 1);
y= @(t) -b/2*(1-t);
dx =@(t) 0*ones(numel(t), 1);
dy =@(t) b/2*ones(numel(t), 1);
ddx =@(t) 0*ones(numel(t), 1);
ddy =@(t) 0*ones(numel(t), 1);
c7 = Interface(x,y,dx,dy,ddx,ddy);

%% c8:
x=@(t) a*ones(numel(t), 1);
y= @(t) b/2*(1+t);
dx =@(t) 0*ones(numel(t), 1);
dy =@(t) b/2*ones(numel(t), 1);
ddx =@(t) 0*ones(numel(t), 1);
ddy =@(t) 0*ones(numel(t), 1);
c8 = Interface(x,y,dx,dy,ddx,ddy);

%% c9:
x=@(t) a/2*(1-t);
y= @(t) 0*ones(numel(t), 1);
dx =@(t) -a/2*ones(numel(t), 1);
dy =@(t) 0*ones(numel(t), 1);
ddx =@(t) 0*ones(numel(t), 1);
ddy =@(t) 0*ones(numel(t), 1);
c9 = Interface(x,y,dx,dy,ddx,ddy, 0);

%% c10:
x=@(t) -a/2*(1+t);
y= @(t) 0*ones(numel(t), 1);
dx =@(t) -a/2*ones(numel(t), 1);
dy =@(t) 0*ones(numel(t), 1);
ddx =@(t) 0*ones(numel(t), 1);
ddy =@(t) 0*ones(numel(t), 1);
c10 = Interface(x,y,dx,dy,ddx,ddy, 0);



%% c11:
x=@(t) 0*ones(numel(t), 1);
y= @(t) b/2*(1-t);
dx =@(t) 0*ones(numel(t), 1);
dy =@(t) -b/2*ones(numel(t), 1);
ddx =@(t) 0*ones(numel(t), 1);
ddy =@(t) 0*ones(numel(t), 1);
c11 = Interface(x,y,dx,dy,ddx,ddy, 0);


%% c12:
x=@(t) 0*ones(numel(t), 1);
y= @(t) -b/2*(1+t);
dx =@(t) 0*ones(numel(t), 1);
dy =@(t) -b/2*ones(numel(t), 1);
ddx =@(t) 0*ones(numel(t), 1);
ddy =@(t) 0*ones(numel(t), 1);
c12 = Interface(x,y,dx,dy,ddx,ddy, 0);


% top left.
%NOT USE id =0!!!!
 Interfz = cell(4,1);
 Interfz{1} =c2;
 Interfz{2} =c3;
 Interfz{3} =c10.SetDirection(-1);
 Interfz{4} =c11.SetDirection(-1);
 k1=4;
 d1 =  Domain(1,0,-1.0*k1,Interfz);


% bottom left.
%NOT USE id =0!!!!
 Interfz = cell(4,1);
 Interfz{1} =c10;
 Interfz{2} =c4;
 Interfz{3} =c5;
 Interfz{4} =c12.SetDirection(-1);
 k1=4;
 d2 =  Domain(2,0,-1.0*k1,Interfz);

% bottom right.
%NOT USE id =0!!!!
 Interfz = cell(4,1);
 Interfz{1} =c9;
 Interfz{2} =c12;
 Interfz{3} =c6;
 Interfz{4} =c7;
 k1=4;
 d3 =  Domain(3,0,-1.0*k1,Interfz);

% top right.
%NOT USE id =0!!!!
 Interfz = cell(4,1);
 Interfz{1} =c1;
 Interfz{2} =c11;
 Interfz{3} =c9.SetDirection(-1);
 Interfz{4} =c8;
 k1=4;
 d4 =  Domain(4,0,-1.0*k1,Interfz);



%Exterior domain
Interfz = cell(8, 1);
Interfz{1} = c2.SetDirection(-1);
Interfz{2} = c1.SetDirection(-1);
Interfz{3} = c8.SetDirection(-1);
Interfz{4} = c7.SetDirection(-1);
Interfz{5} = c6.SetDirection(-1);
Interfz{6} = c5.SetDirection(-1);
Interfz{7} = c4.SetDirection(-1);
Interfz{8} = c3.SetDirection(-1);
d0 = Domain(0,0,-1.0*k1,Interfz);



k_hlf = 1i*p(R*omega.^(-(0:floor(M/2))))/dt;
[xg,wg] = lgwt(Ng,-1,1);
Un =@(n,t) sin((n+1)*acos(t))./sin(acos(t));


%% Evaluate the field

pause(2)
% 
Xlim = [-1 1];
Ylim = [-1 1];

xx = linspace(Xlim(1),Xlim(2),floor(abs(Xlim(2)-Xlim(1))*50));
yy = linspace(Ylim(1),Ylim(2),floor(abs(Ylim(2)-Ylim(1))*50));
[X,Y] = meshgrid(xx,yy);

Nx = numel(xx);
Ny = numel(yy);

% pts = [reshape(X,Nx*Ny,1) reshape(Y,Nx*Ny,1)];
% 
% pts = [1 0;-0.25 0;0.25 0];

pts = [b/2 a/2;-b/2 a/2;-b/2 -a/2;b/2 -a/2; 2*b 2*a];
Npts = size(pts, 1);

dom0 = [b a;-b a;-b -a;b -a;b a];
dom1 = [0 a;-b a;-b 0;0 0;0 a];
dom2 = [0 0;-b 0;-b -a;0 -a;0 0];
dom3 = [b 0;0 0;0 -a;b -a;b 0];
dom4 = [b a;0 a;0 0;b 0;b a];

% in1 = inpolygon(pts(:,1),pts(:,2),vert_1(:,1),vert_1(:,2));
in1 = inpolygon(pts(:,1),pts(:,2),dom1(:,1),dom1(:,2));
in2 = inpolygon(pts(:,1),pts(:,2),dom2(:,1),dom2(:,2));
in3 = inpolygon(pts(:,1),pts(:,2),dom3(:,1),dom3(:,2));
in4 = inpolygon(pts(:,1),pts(:,2),dom4(:,1),dom4(:,2));
in0 = inpolygon(pts(:,1),pts(:,2),dom0(:,1),dom0(:,2));
in0 = ~in0;


k_hlf = 1i*p(R*omega.^(-(0:floor(M/2))))/dt;
[xg,wg] = lgwt(Ng,-1,1);
% xg = linspace(-1, 1, Ng+2);
% xg = xg(2:end-1).';
Un =@(n,t) sin((n+1)*acos(t))./sin(acos(t));



disp('Potentials')

traces0 = cell(1,4);
traces = cell(4,4);
lambda = R;
for dom = 1:4
dira = tD(dom, 1, :, :);
neua = tN(dom, 1, :, :);

dira = reshape(dira, Ng, M+1);
neua = reshape(neua, Ng, M+1);

dirb = tD(dom, 2, :, :);
neub = tN(dom, 2, :, :);

dirb = reshape(dirb, Ng, M+1);
neub = reshape(neub, Ng, M+1);

dirc = tD(dom, 3, :, :);
neuc = tN(dom, 3, :, :);

dirc = reshape(dirc, Ng, M+1);
neuc = reshape(neuc, Ng, M+1);

dird = tD(dom, 4, :, :);
neud = tN(dom, 4, :, :);

dird = reshape(dird, Ng, M+1);
neud = reshape(neud, Ng, M+1);

dir = zeros(4*Ng, M+1);
neu = zeros(4*Ng, M+1);

for n = 1:M+1
   dir(:, n) = [dira(:, n); dirb(:, n);dirc(:, n);dird(:, n)];  
   neu(:, n) = [neua(:, n); neub(:, n);neuc(:, n);neud(:, n)];
end

Lam = repmat(lambda.^(0:M),4*Ng,1);
dir = fft(Lam.*dir,[],2);
neu = fft(Lam.*neu,[],2);

traces{dom, 1} = dir;
traces{dom, 2} = neu;
end





for dom = 1:1
dira = tD0(dom, 1, :, :);
neua = tN0(dom, 1, :, :);

dira = reshape(dira, Ng, M+1);
neua = reshape(neua, Ng, M+1);

dirb = tD0(dom, 2, :, :);
neub = tN0(dom, 2, :, :);

dirb = reshape(dirb, Ng, M+1);
neub = reshape(neub, Ng, M+1);

dirc = tD0(dom, 3, :, :);
neuc = tN0(dom, 3, :, :);

dirc = reshape(dirc, Ng, M+1);
neuc = reshape(neuc, Ng, M+1);

dird = tD0(dom, 4, :, :);
neud = tN0(dom, 4, :, :);

dird = reshape(dird, Ng, M+1);
neud = reshape(neud, Ng, M+1);


dire = tD0(dom, 5, :, :);
neue = tN0(dom, 5, :, :);

dire = reshape(dire, Ng, M+1);
neue = reshape(neue, Ng, M+1);

dirf = tD0(dom, 6, :, :);
neuf = tN0(dom, 6, :, :);

dirf = reshape(dirf, Ng, M+1);
neuf = reshape(neuf, Ng, M+1);

dirg = tD0(dom, 7, :, :);
neug = tN0(dom, 7, :, :);

dirg = reshape(dirg, Ng, M+1);
neug = reshape(neug, Ng, M+1);

dirh = tD0(dom, 8, :, :);
neuh = tN0(dom, 8, :, :);

dirh = reshape(dirh, Ng, M+1);
neuh = reshape(neuh, Ng, M+1);


dir = zeros(8*Ng, M+1);
neu = zeros(8*Ng, M+1);

for n = 1:M+1
   dir(:, n) = [dira(:, n); dirb(:, n);dirc(:, n);dird(:, n); dire(:, n); dirf(:, n);dirg(:, n);dirh(:, n)];  
   neu(:, n) = [neua(:, n); neub(:, n);neuc(:, n);neud(:, n); neue(:, n); neuf(:, n);neug(:, n);neuh(:, n)];
end

Lam = repmat(lambda.^(0:M),8*Ng,1);
dir = fft(Lam.*dir,[],2);
neu = fft(Lam.*neu,[],2);

traces0{dom, 1} = dir;
traces0{dom, 2} = neu;
end

%%

u_hlf = zeros(Npts,M/2+1);

hbar = parfor_progressbar(M/2+1,'Compute potentials...'); %create the progress bar 
tic

c = [c_0 c_1 c_2 c_3 c_4];
for n=1:M/2+1

    k = k_hlf(n)./c;


%     mu = phip_hlf(:, n);

    phi0 = traces0{1, 1}(:, n);
    psi0 = traces0{1, 2}(:, n);
    
    phi1 = traces{1, 1}(:, n);
    psi1 = traces{1, 2}(:, n);
    
    phi2 = traces{2, 1}(:, n);
    psi2 = traces{2, 2}(:, n);
    
    phi3 = traces{3, 1}(:, n);
    psi3 = traces{3, 2}(:, n);
    
    phi4 = traces{4, 1}(:, n);
    psi4 = traces{4, 2}(:, n);
    
    
    u0 = zeros(sum(in0), 1);    
    u1 = zeros(sum(in1), 1);
    u2 = zeros(sum(in2), 1);
    u3 = zeros(sum(in3), 1);
    u4 = zeros(sum(in4), 1);
    

    pot0 = compute_potentials(k(1),pts(in0,:),d0, Ng);

    pot1 = compute_potentials(k(2),pts(in1,:),d1, Ng);

    pot2 = compute_potentials(k(3),pts(in2,:),d2, Ng);
    
    pot3 = compute_potentials(k(4),pts(in3,:),d3, Ng);

    pot4 = compute_potentials(k(5),pts(in4,:),d4, Ng);

    u0 = pot0.SL*psi0-pot0.DL*phi0;

    u1 = pot1.SL*psi1-pot1.DL*phi1;
    
    u2 = pot2.SL*psi2-pot2.DL*phi2;
    
    u3 = pot3.SL*psi3-pot3.DL*phi3;
    
    u4 = pot4.SL*psi4-pot4.DL*phi4;
    % Putting everything together
    u_hlf(in0, n) = u0;
    u_hlf(in1, n) = u1;
    u_hlf(in2, n) = u2; 
    u_hlf(in3, n) = u3;
    u_hlf(in4, n) = u4; 
    hbar.iterate(1);
toc
end
toc
close(hbar);

% up = [u_hlf(:,1) conj(u_hlf(:,2:end)) u_hlf(:,end:-1:2)];


u_hlf(:,M+2-(1:floor(M/2)))=conj(u_hlf(:,2:floor(M/2)+1));  % mirror the hermitian seq
u=real(ifft(u_hlf,[],2));                            % w
u=bsxfun(@times,u,R.^(-(0:M)));                % g

u0 = u(in0, :);
u1 = u(in1, :);
u2 = u(in2, :);
u3 = u(in3, :);
u4 = u(in4, :);

clear u dir neu traces



if eme == 0
    u0_Ref = u0;
    u1_Ref = u1;
    u2_Ref = u2;
    u3_Ref = u3;
    u4_Ref = u4;
    tD_Ref = tD;
    tN_Ref = tN;
    tD0_Ref = tD0;
    tN0_Ref = tN0;
else
disp('Errors')


load 'Ai_square.mat'
cte = 101*8;

A0 = Ai(1:2*cte, 1:2*cte); A0 = A0(1:size(A0,1)/2, size(A0, 1)/2+1:end);
A1 = Ai(2*cte + (1:cte),2*cte + (1:cte)); A1 = A1(1:size(A1,1)/2, size(A1, 1)/2+1:end);
A2 = Ai(3*cte + (1:cte),3*cte + (1:cte)); A2 = A2(1:size(A2,1)/2, size(A2, 1)/2+1:end);
A3 = Ai(4*cte + (1:cte),4*cte + (1:cte)); A3 = A3(1:size(A3,1)/2, size(A3, 1)/2+1:end);
A4 = Ai(5*cte + (1:cte),5*cte + (1:cte)); A4 = A4(1:size(A4,1)/2, size(A4, 1)/2+1:end);

errorD(eme, 1) = traceErrorDom_square(A0, tD0_Ref, tD0,0, 1);
errorD(eme, 2) = traceErrorDom_square(A1, tD_Ref, tD,0, 2);
errorD(eme, 3) = traceErrorDom_square(A2, tD_Ref, tD,0, 3);
errorD(eme, 4) = traceErrorDom_square(A3, tD_Ref, tD,0, 4);
errorD(eme, 5) = traceErrorDom_square(A4, tD_Ref, tD,0, 5)


errorN(eme, 1) = traceErrorDom_square(A0, tN0_Ref, tN0,1, 1);
errorN(eme, 2) = traceErrorDom_square(A1, tN_Ref, tN,1, 2);
errorN(eme, 3) = traceErrorDom_square(A2, tN_Ref, tN,1, 3);
errorN(eme, 4) = traceErrorDom_square(A3, tN_Ref, tN,1, 4);
errorN(eme, 5) = traceErrorDom_square(A4, tN_Ref, tN,1, 5)

ind = 1:Mmax/M:Mmax+1;
errorU(eme, 1) = norm(u0_Ref(ind) - u0)/norm(u0_Ref(ind)); 
errorU(eme, 2) = norm(u1_Ref(ind) - u1)/norm(u1_Ref(ind)); 
errorU(eme, 3) = norm(u2_Ref(ind) - u2)/norm(u2_Ref(ind)); 
errorU(eme, 4) = norm(u3_Ref(ind) - u3)/norm(u3_Ref(ind));
errorU(eme, 5) = norm(u4_Ref(ind) - u4)/norm(u4_Ref(ind))


end
clearvars inc inc2;


eme = eme +1;

save(['square-solution-',num2str(M),'.mat']);




end

% % eNes = eNes(2:end);
% subplot(2, 1, 1);
% eNes = emes;
% loglog(eNes, errorD(:,1), '-x', ...
%     eNes, errorD(:,2), '-o', ...
%     eNes, errorD(:,3), '-s', ...
%     eNes, errorD(:,4), '-^', ...
%     eNes, errorD(:,5), '-d', ...
%     eNes, errorD(:,6), '-+', ...
%     eNes, 1e5*eNes.^(-2), '--');
% %%
% subplot(2, 1, 2);
% loglog(...
%     eNes, errorN(:,1), '-x', ...
%     eNes, errorN(:,2), '-o', ...
%     eNes, errorN(:,3), '-s', ...
%     eNes, errorN(:,4), '-^', ...
%     eNes, errorN(:,5), '-d', ...
%     eNes, errorN(:,6), '-+', ...
%     eNes, 1.8e2*eNes.^(-1), '--');



save('test.mat', 'errorD', 'errorN', 'errorU');
% mail_results('Resultados TDMTF', date);
%%
% emes = [50 100 200 400];
% % % subplot(2, 1, 1);
loglog(emes, errorD(:, 1), '-s');hold on;
loglog(emes, errorD(:, 2), '-^');
loglog(emes, errorD(:, 3), '-d');
loglog(emes, errorD(:, 4), '-^');
loglog(emes, errorD(:, 5), '-d');
loglog(emes, 11*emes.^(-3/2), '--k');
loglog(emes, 0.75e4*emes.^(-2), '-.k');
legend({'dom0', 'dom1', 'dom2','Order 3/2', 'Order 2'},...
        'Location', 'southwest');hold off;
% 


saveas(gcf, 'errorD.png');
close all
% % % % % subplot(2, 1, 2);    
loglog(emes, errorN(:, 1), '-s');hold on;
loglog(emes, errorN(:, 2), '-^');
loglog(emes, errorN(:, 3), '-d');
loglog(emes, errorN(:, 4), '-^');
loglog(emes, errorN(:, 5), '-d');
loglog(emes, 55*emes.^(-1), '--k');
loglog(emes, 0.75e2*emes.^(-3/2), '-.k');
legend({'dom0', 'dom1', 'dom2','Order 1/2', 'Order 1'},...
        'Location', 'southwest');

saveas(gcf, 'errorN.png');
close all
% 
% 
% % 
loglog(emes, errorU(:, 1), '-s');hold on;
loglog(emes, errorU(:, 2), '-^');
loglog(emes, errorU(:, 3), '-d');
loglog(emes, errorU(:, 4), '-^');
loglog(emes, errorU(:, 5), '-d');
loglog(emes, 11*emes.^(-3/2), '--k');
loglog(emes, 0.75e4*emes.^(-2), '-.k');
legend({'0', '1', '2','Order 1', 'Order 2'},...
        'Location', 'southwest');hold off;


saveas(gcf, 'errorU.png');
close all

mail_results('TDMTF', 'Results', {'errorD.png', 'errorN.png', 'errorU.png'})