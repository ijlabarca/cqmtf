

clear
close all
clc

%%
%-----------------------------------------
% set main parameters
%-----------------------------------------
eme = 0;
ene = 0;


M = 200; % timesteps
T = 10;
dt = T/M;
time = 0:dt:T;
% time = linspace(0, T, M);
tlag = 0.5;
dirx = sqrt(0.5);
diry = -sqrt(0.5);
disp('-----------------------------------------------')
% wavespeeds
Nmax = 14;
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
        


%     if norm(solver.m_b(:,l+1)) > 1e-12
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

    g(:,l+1)= solver.m_A \ solver.m_b(:,l+1);  
    
%     sol = solver.m_A \ solver.m_b(:,l+1);  

%     save(['solutions/solution-', num2str(l), '.mat'], 'sol');
    
%     clear sol 
%     disp(solver.CalderonError(g(:, l+1)) );


%     if l == 0
%        toc 
%     end

%     [g(:,l+1), ~] = gmres(solver.m_A, solver.m_b(:,l+1));
% 
%             disp(k0);
%     solver.JumpError(g(:, l+1), l, c_0, T, M, tlag, dirx, diry, d(1))

%     norm(g(:,l+1))

%     end
    
    hbar.iterate(1);
toc
end
toc
close(hbar);
disp(['Number of timesteps: ', num2str(M)]);

disp(' ')
% 
% for l=0:M/2
%     load(['solutions/solution-', num2str(l), '.mat']);
%     g(:, l+1) = sol;
%     
%     clear sol 
%     
%     
% end


g(:,M+2-(1:floor(M/2)))=conj(g(:,2:floor(M/2)+1));  % mirror the hermitian seq
g=real(ifft(g,[],2));                            % w
g=bsxfun(@times,g,R.^(-(0:M)));                % g




%%

disp('Computing Traces');
% Ng = 2*max(solver.m_N{1});

Ng = 101;
[xg,wg] = lgwt(Ng,-1,1);
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

lambdaD0(1, :, :) = g(1:2*dim, :); g(1:2*dim, :) = [];
lambdaN0(1, :, :) = g(1:2*dim, :); g(1:2*dim, :) = [];


lambdaD(1, :, :) = g(1:dim, :); g(1:dim, :) = [];
lambdaN(1, :, :) = g(1:dim, :); g(1:dim, :) = [];


lambdaD(2, :, :) = g(1:dim, :); g(1:dim, :) = [];
lambdaN(2, :, :) = g(1:dim, :); g(1:dim, :) = [];

lambdaD(3, :, :) = g(1:dim, :); g(1:dim, :) = [];
lambdaN(3, :, :) = g(1:dim, :); g(1:dim, :) = [];

lambdaD(4, :, :) = g(1:dim, :); g(1:dim, :) = [];
lambdaN(4, :, :) = g(1:dim, :); g(1:dim, :) = [];
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



%% % Plot Traces
% % 
% plotTraces(deltat, tD,0,uinc, 1, 2, 1, 1);
% 
% 
% plotTraces(deltat, tD,0,uinc, 2, 3, 2, 1);
% 
% 
% plotTraces(deltat, tD,0,uinc, 1, 3, 2, 2);



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


%  % We save all the domains in one cell array
%  Domains = cell(4,1);
%  Domains{1} =d1;
%  Domains{2} =d2;
%  Domains{3} =d3;
%  Domains{4} =d4;


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

% % We save all the domains in one cell array
% Domains = cell(2,1);
% Domains{1} =d1;
% Domains{2} =d2;

%     M = 100;


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

pts = [reshape(X,Nx*Ny,1) reshape(Y,Nx*Ny,1)];



% % vert_1 = [x1(end:-1:1,:);[x1(1,1) h2];[x1(end,1) 100]];
% [x1, y1] = d0.geo(xg, 1);
% [x2, y2] = d0.geo(xg, 2);
% [x3, y3] = d0.geo(xg, 3);
% [x4, y4] = d0.geo(xg, 4);
% [x5, y5] = d0.geo(xg, 5);
% [x6, y6] = d0.geo(xg, 6);
% [x7, y7] = d0.geo(xg, 7);
% [x8, y8] = d0.geo(xg, 8);
% dom0 = [x1 y1;x2 y2;x3 y3;x4 y4;x5 y5;x6 y6;x7 y7;x8 y8];

dom0 = [b a;-b a;-b -a;b -a;b a];


% [x1, y1] = d1.geo(xg, 1);
% [x2, y2] = d1.geo(xg, 2);
% [x3, y3] = d1.geo(xg, 3);
% [x4, y4] = d1.geo(xg, 4);
% dom1 = [x1 y1;x2 y2;x3 y3;x4 y4];

dom1 = [0 a;-b a;-b 0;0 0;0 a];

% [x1, y1] = d2.geo(xg, 1);
% [x2, y2] = d2.geo(xg, 2);
% [x3, y3] = d2.geo(xg, 3);
% [x4, y4] = d2.geo(xg, 4);
% dom2 = [x1 y1;x2 y2;x3 y3;x4 y4];

dom2 = [0 0;-b 0;-b -a;0 -a;0 0];

% [x1, y1] = d3.geo(xg, 1);
% [x2, y2] = d3.geo(xg, 2);
% [x3, y3] = d3.geo(xg, 3);
% [x4, y4] = d3.geo(xg, 4);
% dom3 = [x1 y1;x2 y2;x3 y3;x4 y4];

dom3 = [b 0;0 0;0 -a;b -a;b 0];

% [x1, y1] = d4.geo(xg, 1);
% [x2, y2] = d4.geo(xg, 2);
% [x3, y3] = d4.geo(xg, 3);
% [x4, y4] = d4.geo(xg, 4);
% dom4 = [x1 y1;x2 y2;x3 y3;x4 y4];

dom4 = [b a;0 a;0 0;b 0;b a];



% in1 = inpolygon(pts(:,1),pts(:,2),vert_1(:,1),vert_1(:,2));
in1 = inpolygon(pts(:,1),pts(:,2),dom1(:,1),dom1(:,2));
in2 = inpolygon(pts(:,1),pts(:,2),dom2(:,1),dom2(:,2));
in3 = inpolygon(pts(:,1),pts(:,2),dom3(:,1),dom3(:,2));
in4 = inpolygon(pts(:,1),pts(:,2),dom4(:,1),dom4(:,2));
in0 = inpolygon(pts(:,1),pts(:,2),dom0(:,1),dom0(:,2));
in0 = ~in0;
%%
% Incident field



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

u_hlf = zeros(Nx*Ny,M/2+1);

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


% Inverting Z-transform
% Lam = repmat(lambda.^(0:M),size(pts,1),1);
% 
% u = Lam.^(-1).*ifft(up,[],2);

uinc = zeros(sum(in0), M+1);
z = pts(in0, 1) + 1i*pts(in0, 2);
for m = 1:M+1
onda = timedomain_wavefunction2d((m-1)*deltat, tlag, dirx, diry, c_0);
uinc(:, m) = onda.evaluate(z);

clear onda;
                    
end


u(in0,:) =  u(in0,:) + uinc;    


%% Plot the solution
close all  

X = reshape(pts(:,1),Ny,Nx);

Y = reshape(pts(:,2),Ny,Nx);
% % 
% vid =  VideoWriter('composite.avi');

% open(vid)
figure;

count = 1;
for n=1:M+1
    

    if mod(n, 1) == 0

%         disp(n);
%     contourf(xx,yy,reshape(real(u(:, n)),Ny,Nx), 20); hold on;
    surf(X,Y,reshape(real(u(:, n)),Ny,Nx),'EdgeColor','none'); hold on;
    plot3(dom0(:,1),dom0(:,2),10*ones(size(dom0,1),1),'LineWidth',4,'color', 'k');
    plot3(dom1(:,1),dom1(:,2),10*ones(size(dom1,1),1),'LineWidth',4,'color', 'k');
    plot3(dom2(:,1),dom2(:,2),10*ones(size(dom2,1),1),'LineWidth',4,'color', 'k');
    plot3(dom3(:,1),dom3(:,2),10*ones(size(dom3,1),1),'LineWidth',4,'color', 'k');
    plot3(dom4(:,1),dom4(:,2),10*ones(size(dom4,1),1),'LineWidth',4,'color', 'k');
%     title(['Plane wave, ','$\theta = \pi/4$'],'Interpreter','latex');
% title((n-1)*dt);
    view(2);
    colormap(brewermap([],'*RdBu'))
%     colormap(brewermap([],'*RdGy'))
    caxis([-0.7 0.7])
    shading interp;
    hold off;
%     colorbar
    axis equal
    axis([Xlim Ylim])
% %     
%     frm = getframe(gcf);
%     writeVideo(vid,frm)
%     

    set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.
    set(gca,'Yticklabel',[]) %to just get rid of the numbers but leave the ticks.

drawnow;


% 

    
    count = count+1;
    pause(0.1)
%     disp(deltat*(n-1));
    end

end