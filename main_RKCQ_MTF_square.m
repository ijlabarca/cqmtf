
%% RUNGE KUTTA CONVOLUTION QUADRATURE

clear 
close all
clc

%% Runge-Kutta parameters
Ark=[5/12 -1/12; 3/4 1/4];   % RADAU IIa
% Ark=[1/6 -1/3 1/6; 1/6 5/12 -1/12; 1/6 2/3 1/6]; % LOBATTO IIIc
bb=Ark(end,:);
S=size(Ark,1); RK = S;
Am1 = inv(Ark);
crk=Ark*ones(S,1);
B = (Am1*ones(size(Ark,1),1))*[zeros(1,S-1),1];

%% PDE parameters
M = 800; % timesteps
disp('-----------------------------------------------')
% wavespeeds
Nmax = 10;
% M = emes(eme);


c_0 = 1.;
c_1 = 0.5;
c_2 = 0.25;
c_3 = 0.5;
c_4 = 0.25;


%% Domain
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
 
dim = 24*2*(2*Nmax + 1);
idx=@(s) (s-1)*dim+1:s*dim;



T = 20;
dt = T/M;
time = 0:dt:T;
tlag = 0.5;
dirx = sqrt(0.5);
diry = -sqrt(0.5);

omega = exp(2*pi*1i/(M));
R = eps^(0.5/(M));

%% Right hand side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% inc{1} = tdwavefunction2d(T-dt, M-1, 1,tlag,dirx, diry, c_0, crk(1)*dt, 1);
% 
% inc2{1} = tdwavefunction2d(T-dt, M-1, 1,tlag,dirx, diry, c_0, crk(2)*dt, 1);
% 
% inc3{1} = tdwavefunction2d(T-dt, M-1, 1,tlag,dirx, diry, c_0, crk(3)*dt, 1);


inc{1} = td_test(T-dt, M-1, 1,tlag,dirx, diry, c_0,1, crk(1)*dt, 1);

inc2{1} = td_test(T-dt, M-1, 1,tlag,dirx, diry, c_0,1, crk(2)*dt, 1);

% inc3{1} = td_test(T-dt, M-1, 1,tlag,dirx, diry, c_0,1, crk(3)*dt, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup solver with empty incident field for now
kwave_real = 0;
kwave_imag = -4*pi;


% Nmax = ceil(1.1*maxK)+7
%%

solver = solverMtf28(c_0,kwave_real,kwave_imag,inc,Domains);
solver2 = solverMtf28(c_0,kwave_real,kwave_imag,inc2,Domains);
% solver3 = solverMtf28(c_0,kwave_real,kwave_imag,inc3,Domains);

disp('Lado derecho')
tic
solver.setupb(Nmax, M-1)
toc
solver2.setupb(Nmax, M-1)
toc
% solver3.setupb(Nmax, M-1)
toc

F =zeros(dim,S,M);
F(:, 1, :) = solver.m_b;
F(:, 2, :) = solver2.m_b;
% F(:, 3, :) = solver3.m_b;
F = reshape(F, [S*dim, M]);

g = zeros(dim*S,M);

% load 'sol-ram.mat'
%%
disp(' ')
disp('TDMTF-solve');
tic
hbar = parfor_progressbar(M,'Solving Linear Systems...'); %create the progress bar 
for l=0:M-1

    [P,Lambda]=eig(Am1-R*omega^(-(l))*B);   
    Lambda=diag(Lambda)/dt;
    gl=kron(inv(P),sparse(eye(dim)))*F(:,l+1);

    ul=zeros(S*dim,1);    
    
    for s = 1:S
    k_aux = Lambda(s)./[c_0; c_1; c_2; c_3; c_4];
    k0 = k_aux(1);
    k1 = k_aux(2);
    k2 = k_aux(3);
    k3 = k_aux(4);
    k4 = k_aux(5);

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
    
    rhs = gl(idx(s));
    
    ul(idx(s))= solver.m_A \ rhs;  
    end
    
    g(:, l+1) = kron(P,sparse(eye(dim)))*ul;
    
    hbar.iterate(1);
    toc
end
toc
close(hbar);
disp(['Number of timesteps: ', num2str(M)]);

disp(' ')

g=real(ifft(g,[],2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%

g=bsxfun(@times,g,R.^(-(0:M-1)));
disp('Computing Traces');


%%
Ng = 101;
[xg,wg] = lgwt(Ng,-1,1);
Un =@(n,t) sin((n+1)*acos(t))./sin(acos(t));


tD0 = zeros(1*S,8, Ng, M);
tN0 = zeros(1*S,8, Ng, M);


tD = zeros(4*S,4, Ng, M);
tN = zeros(4*S,4, Ng, M);

deltat = T/M;

uinc = zeros(8,Ng, M);
dnuinc = zeros(8,Ng, M);

dim_intf = 4*(2*Nmax +1);
lambdaD0 = zeros(S, 2*dim_intf, M);
lambdaN0 = zeros(S, 2*dim_intf, M);

lambdaD = zeros(S*4, dim_intf, M);
lambdaN = zeros(S*4, dim_intf, M);

for s = 1:S
gg = g(idx(s), :);

lambdaD0(s, :, :) = gg(1:2*dim_intf, :); gg(1:2*dim_intf, :)=[];
lambdaN0(s, :, :) = gg(1:2*dim_intf, :); gg(1:2*dim_intf, :)=[];

lambdaD((s-1)*4 + 1, :, :) = gg(1:dim_intf, :); gg(1:dim_intf, :)=[];
lambdaN((s-1)*4 + 1, :, :) = gg(1:dim_intf, :); gg(1:dim_intf, :)=[];

lambdaD((s-1)*4 + 2, :, :) = gg(1:dim_intf, :); gg(1:dim_intf, :)=[];
lambdaN((s-1)*4 + 2, :, :) = gg(1:dim_intf, :); gg(1:dim_intf, :)=[];

lambdaD((s-1)*4 + 3, :, :) = gg(1:dim_intf, :); gg(1:dim_intf, :)=[];
lambdaN((s-1)*4 + 3, :, :) = gg(1:dim_intf, :); gg(1:dim_intf, :)=[];

lambdaD((s-1)*4 + 4, :, :) = gg(1:dim_intf, :); gg(1:dim_intf, :)=[];
lambdaN((s-1)*4 + 4, :, :) = gg(1:dim_intf, :); gg(1:dim_intf, :)=[];

clear gg;
end
%%
tic




for s = 1:S
for dom = 1:5


    counter = 0;

     NIntrfz0 = solver.m_DomainArray{dom}.m_NumInfefaces;

     for jj=1:NIntrfz0



               [nx,ny] = solver.m_DomainArray{dom}.normal(xg,jj);

               [x,y] = solver.m_DomainArray{dom}.geo(xg,jj);

               z = x + 1i*y;

               J = solver.m_DomainArray{dom}.J(xg, jj);



               for m = 1:M


                   for n = 0:2*solver.m_N{dom}(jj)


                        for gg=1:Ng



                          if dom == 1
                           if length(J) == 1


                               tD0(s,jj, gg, m) = tD0(s,jj, gg, m) + (lambdaD0(s, counter + n + 1, m) * Un(n, xg(gg))./J);  

                               tN0(s,jj, gg, m) = tN0(s,jj, gg, m) + (lambdaN0(s, counter + n + 1, m) * Un(n, xg(gg))./J); 

                           else 

                               tD0(s,jj, gg, m) = tD0(s,jj, gg, m) + (lambdaD0(s, counter + n + 1, m) * Un(n, xg(gg))./J(gg));  

                               tN0(s,jj, gg, m) = tN0(s,jj, gg, m) + (lambdaN0(s, counter + n + 1, m) * Un(n, xg(gg))./J(gg)); 


                           end
                          else
                              
                            if length(J) == 1


                               tD((s-1)*4 + dom-1,jj, gg, m) = tD((s-1)*4 + dom-1,jj, gg, m) + (lambdaD((s-1)*4 + dom-1, counter + n + 1, m) * Un(n, xg(gg))./J);  

                               tN((s-1)*4 + dom-1,jj, gg, m) = tN((s-1)*4 + dom-1,jj, gg, m) + (lambdaN((s-1)*4 + dom-1, counter + n + 1, m) * Un(n, xg(gg))./J); 

                            else 

                               tD((s-1)*4 + dom-1,jj, gg, m) = tD((s-1)*4 + dom-1,jj, gg, m) + (lambdaD((s-1)*4 + dom-1, counter + n + 1, m) * Un(n, xg(gg))./J(gg));  

                               tN((s-1)*4 + dom-1,jj, gg, m) = tN((s-1)*4 + dom-1,jj, gg, m) + (lambdaN((s-1)*4 + dom-1, counter + n + 1, m) * Un(n, xg(gg))./J(gg)); 


                            end
                           
                                                         
                          end
                           
                           
                           
                           
                           
                        end


                   end


                   if dom == 1 && s == S
                       onda = timedomain_wavefunction2d(m*deltat, tlag, dirx, diry, c_0);
                       [dx, dy] = onda.evaluateGradient(z);

                       dnuinc(jj,:,m) = nx .* dx + ny .* dy;
                       uinc(jj,:, m) = onda.evaluate(z);

                   end


               end

               counter = counter + 2*solver.m_N{dom}(jj) + 1;

     end

end

end

toc



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





[xg,wg] = lgwt(Ng,-1,1);
Un =@(n,t) sin((n+1)*acos(t))./sin(acos(t));


%% Evaluate the field

pause(2)
% 
Xlim = [-1 1];
Ylim = [-1 1];

xx = linspace(Xlim(1),Xlim(2),floor(abs(Xlim(2)-Xlim(1))*80));
yy = linspace(Ylim(1),Ylim(2),floor(abs(Ylim(2)-Ylim(1))*80));
[X,Y] = meshgrid(xx,yy);

Nx = numel(xx);
Ny = numel(yy);

pts = [reshape(X,Nx*Ny,1) reshape(Y,Nx*Ny,1)];


Npts = size(pts, 1);




dom0 = [b a;-b a;-b -a;b -a;b a];
dom1 = [0 a;-b a;-b 0;0 0;0 a];
dom2 = [0 0;-b 0;-b -a;0 -a;0 0];
dom3 = [b 0;0 0;0 -a;b -a;b 0];
dom4 = [b a;0 a;0 0;b 0;b a];


in1 = inpolygon(pts(:,1),pts(:,2),dom1(:,1),dom1(:,2));
in2 = inpolygon(pts(:,1),pts(:,2),dom2(:,1),dom2(:,2));
in3 = inpolygon(pts(:,1),pts(:,2),dom3(:,1),dom3(:,2));
in4 = inpolygon(pts(:,1),pts(:,2),dom4(:,1),dom4(:,2));
in0 = inpolygon(pts(:,1),pts(:,2),dom0(:,1),dom0(:,2));
in0 = ~in0;


%% POTENTIALS

disp('Potentials')
traces0 = cell(S,8);
traces = cell(4*S,4);

lambda = R;
for s = 1:S
for dom = 1:4
dira = tD((s-1)*4 + dom, 1, :, :);
neua = tN((s-1)*4 + dom, 1, :, :);

dira = reshape(dira, Ng, M);
neua = reshape(neua, Ng, M);

dirb = tD((s-1)*4 + dom, 2, :, :);
neub = tN((s-1)*4 + dom, 2, :, :);

dirb = reshape(dirb, Ng, M);
neub = reshape(neub, Ng, M);

dirc = tD((s-1)*4 + dom, 3, :, :);
neuc = tN((s-1)*4 + dom, 3, :, :);

dirc = reshape(dirc, Ng, M);
neuc = reshape(neuc, Ng, M);

dird = tD((s-1)*4 + dom, 4, :, :);
neud = tN((s-1)*4 + dom, 4, :, :);

dird = reshape(dird, Ng, M);
neud = reshape(neud, Ng, M);

dir = zeros(4*Ng, M);
neu = zeros(4*Ng, M);

for n = 1:M
   dir(:, n) = [dira(:, n); dirb(:, n); dirc(:, n); dird(:, n)];  
   neu(:, n) = [neua(:, n); neub(:, n); neuc(:, n); neud(:, n)];
end

Lam = repmat(R.^(0:M-1),4*Ng,1);
dir = fft(Lam.*dir,[],2);
neu = fft(Lam.*neu,[],2);

traces{(s-1)*4 + dom, 1} = dir;
traces{(s-1)*4 + dom, 2} = neu;
end

end


for s = 1:S
for dom = 1:1
dira = tD0(s, 1, :, :);
neua = tN0(s, 1, :, :);

dira = reshape(dira, Ng, M);
neua = reshape(neua, Ng, M);

dirb = tD0(s, 2, :, :);
neub = tN0(s, 2, :, :);

dirb = reshape(dirb, Ng, M);
neub = reshape(neub, Ng, M);

dirc = tD0(s, 3, :, :);
neuc = tN0(s, 3, :, :);

dirc = reshape(dirc, Ng, M);
neuc = reshape(neuc, Ng, M);

dird = tD0(s, 4, :, :);
neud = tN0(s, 4, :, :);

dird = reshape(dird, Ng, M);
neud = reshape(neud, Ng, M);


dire = tD0(s, 5, :, :);
neue = tN0(s, 5, :, :);

dire = reshape(dire, Ng, M);
neue = reshape(neue, Ng, M);

dirf = tD0(s, 6, :, :);
neuf = tN0(s, 6, :, :);

dirf = reshape(dirf, Ng, M);
neuf = reshape(neuf, Ng, M);

dirg = tD0(s, 7, :, :);
neug = tN0(s, 7, :, :);

dirg = reshape(dirg, Ng, M);
neug = reshape(neug, Ng, M);

dirh = tD0(s, 8, :, :);
neuh = tN0(s, 8, :, :);

dirh = reshape(dirh, Ng, M);
neuh = reshape(neuh, Ng, M);


dir = zeros(8*Ng, M);
neu = zeros(8*Ng, M);

for n = 1:M
   dir(:, n) = [dira(:, n); dirb(:, n);dirc(:, n);dird(:, n); dire(:, n); dirf(:, n);dirg(:, n);dirh(:, n)];  
   neu(:, n) = [neua(:, n); neub(:, n);neuc(:, n);neud(:, n); neue(:, n); neuf(:, n);neug(:, n);neuh(:, n)];
end

Lam = repmat(lambda.^(0:M-1),8*Ng,1);
dir = fft(Lam.*dir,[],2);
neu = fft(Lam.*neu,[],2);

traces0{s, 1} = dir;
traces0{s, 2} = neu;
end
end







ttraces = [];
for s = 1:S
    aux = [traces0{s, 1}; traces0{s, 2}; ...
           traces{(s-1)*4 + 1, 1};traces{(s-1)*4 + 1, 2}; ...
           traces{(s-1)*4 + 2, 1};traces{(s-1)*4 + 2, 2}; ...
           traces{(s-1)*4 + 3, 1};traces{(s-1)*4 + 3, 2}; ...
           traces{(s-1)*4 + 4, 1};traces{(s-1)*4 + 4, 2}];
    ttraces = [ttraces; aux];
    clear aux;
end



%%
u_hlf = zeros(S*Npts,M);

idx=@(s) (s-1)*Npts+1:(s*Npts);
dim = 24*2*Ng;
idy=@(s) (s-1)*dim+1:(s*dim); 
hbar = parfor_progressbar(M,'Compute potentials...'); %create the progress bar 
tic

c = [c_0 c_1 c_2 c_3 c_4];
for n=1:M
    
    [P,Lambda]=eig(Am1-R*omega^(-(n-1))*B);
    Lambda=diag(Lambda)/dt;    
    
    gl=kron(inv(P),speye(dim))*ttraces(:,n);
      
    ul=zeros(S*Npts,1);
    
    
    for s = 1:S

    u_s = zeros(Npts, 1);
    
    
    k = 1i*Lambda(s)./c;


    mu = gl(idy(s));
    
    phi0 = mu(1:8*Ng); mu(1:8*Ng) = [];
    psi0 = mu(1:8*Ng); mu(1:8*Ng) = [];
    
    phi1 = mu(1:4*Ng); mu(1:4*Ng) = [];
    psi1 = mu(1:4*Ng); mu(1:4*Ng) = [];
    
    phi2 = mu(1:4*Ng); mu(1:4*Ng) = [];
    psi2 = mu(1:4*Ng); mu(1:4*Ng) = [];
    
    phi3 = mu(1:4*Ng); mu(1:4*Ng) = [];
    psi3 = mu(1:4*Ng); mu(1:4*Ng) = [];
    
    phi4 = mu(1:4*Ng); mu(1:4*Ng) = [];
    psi4 = mu(1:4*Ng); mu(1:4*Ng) = [];
    
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
    u_s(in0) = u0;
    u_s(in1) = u1;
    u_s(in2) = u2; 
    u_s(in3) = u3;
    u_s(in4) = u4;
    
    
    ul(idx(s))=u_s;
    end
    
    
    u_hlf(:,n)=kron(P,speye(Npts))*ul;
    
    hbar.iterate(1);
toc
end
toc
close(hbar);

Lam = repmat(R.^(0:M-1),S*Npts,1);

u = Lam.^(-1).*ifft(u_hlf,[],2);

u = u((S-1)*Npts+1:end,:);

uinc = zeros(sum(in0), M);
z = pts(in0, 1) + 1i*pts(in0, 2);
for m = 1:M
onda = timedomain_wavefunction2d(m*deltat, tlag, dirx, diry, c_0);
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
for n=1:M
    

    if mod(n, 1) == 0

%         disp(n);
%     contourf(xx,yy,reshape(real(u(:, n)),Ny,Nx), 20); hold on;
    surf(X,Y,reshape(real(u(:, n)),Ny,Nx),'EdgeColor','none'); hold on;
    plot3(dom0(:,1),dom0(:,2),10*ones(size(dom0,1),1),'LineWidth',2,'color', 'k');
    plot3(dom1(:,1),dom1(:,2),10*ones(size(dom1,1),1),'LineWidth',2,'color', 'k');
    plot3(dom2(:,1),dom2(:,2),10*ones(size(dom2,1),1),'LineWidth',2,'color', 'k');
    plot3(dom3(:,1),dom3(:,2),10*ones(size(dom3,1),1),'LineWidth',2,'color', 'k');
    plot3(dom4(:,1),dom4(:,2),10*ones(size(dom4,1),1),'LineWidth',2,'color', 'k');
    view(2);
    colormap jet
%     colormap(brewermap([],'*RdBu'))
%     colormap(brewermap([],'*RdGy'))
    caxis([-0.8 0.8])
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

    end

end


