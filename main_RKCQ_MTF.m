
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
M = 100; % timesteps
disp('-----------------------------------------------')
% wavespeeds
Nmax = 12;
% M = emes(eme);

c_0 = 1.;
c_1 = 0.5;
c_2 = 0.25;



%% DOMAINs
% origin
center = 0;

% domain wiht one elipse.
a=0.5;
b=0.5;



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

 
% 
% upper half of the kite;
% x = @(t) a*(cos(pi/2*(t+1))+0.65*cos(2*(pi/2*(t+1)))-0.65);
% y = @(t) a*(1.5*sin(pi/2*(t+1)));
% dx = @(t) a*pi/2*(-sin(pi/2*(t+1))-1.3*sin(2*(pi/2*(t+1))));
% dy = @(t) a*pi/2*1.5*cos(pi/2*(t+1));
% ddx = @(t) a*pi/2*pi/2*(-cos(pi/2*(t+1))-2.6*cos(2*(pi/2*(t+1))));
% ddy = @(t) a*pi/2*pi/2*(-1.5*sin(pi/2*(t+1)));
% 
% construction of the curve
% c2 = Interface(x,y,dx,dy,ddx,ddy);
% 
% lower half of the kite;
% x = @(t) a*(cos(pi/2*(t+1)+pi)+0.65*cos(2*(pi/2*(t+1)+pi))-0.65);
% y = @(t) a*(1.5*sin(pi/2*(t+1)+pi));
% dx = @(t) a*pi/2*(-sin(pi/2*(t+1)+pi)-1.3*sin(2*(pi/2*(t+1)+pi)));
% dy = @(t) a*pi/2*1.5*cos(pi/2*(t+1)+pi);
% ddx = @(t) a*pi/2*pi/2*(-cos(pi/2*(t+1)+pi)-2.6*cos(2*(pi/2*(t+1)+pi)));
% ddy = @(t) a*pi/2*pi/2*(-1.5*sin(pi/2*(t+1)+pi));
% c1 = Interface(x,y,dx,dy,ddx,ddy);
% 
% separation between the two halves:
% x=@(t) -a*t;
% y= @(t) 0;
% dx =@(t) -a;
% dy =@(t) 0;
% ddx =@(t) 0;
% ddy =@(t) 0;
% we add a 0 at the end to indicate that this curve is not on the exterior
% domain
% c3 =Interface(x,y,dx,dy,ddx,ddy,0);
% 
% separtion oriented in oposite direction.
% c4 = c3.SetDirection(-1);
% 
% 
% first domain Left domain.
% NOT USE id =0!!!!
%  Interfz = cell(2,1);
%  Interfz{1} =c2;
%  Interfz{2} =c4;
%  k1=4;
%  d1 =  Domain(1,0,-1.0*k1,Interfz);
% 
% 
% 
%  Second domain Right domain.
%  Interfz = cell(2,1);
%  Interfz{1} =c3;
%  Interfz{2} =c1;
%  k2 =8;
%  d2 =  Domain(2,0,-1.0*k2,Interfz);
% 
%  We save all the domains in one cell array
%  Domains = cell(2,1);
%  Domains{1} =d1;
%  Domains{2} =d2;

%%
 
dim = 6*2*(2*Nmax + 1);
idx=@(s) (s-1)*dim+1:s*dim;

T = 20;
dt = T/M;
time = 0:dt:T;
% time = linspace(0, T, M);
tlag = 0.5;
dirx = sqrt(0.5);
diry = -sqrt(0.5);

% dirx = 1;
% diry= 0;
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
    k_aux = Lambda(s)./[c_0; c_1; c_2];
    k0 = k_aux(1);
    k1 = k_aux(2);
    k2 = k_aux(3);

    %Set wave numbers for each subdomain
    solver.m_kext_real = real(k0);

    solver.m_kext_imag = imag(k0);

    solver.m_kext = k0;


    solver.m_DomainArray{1}.m_KwaveNumber = k0;
    solver.m_DomainArray{2}.m_KwaveNumber = k1;
    solver.m_DomainArray{3}.m_KwaveNumber = k2;
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

% xg = linspace(-1, 1, Ng+2);
% xg = xg(2:end-1).';
Un =@(n,t) sin((n+1)*acos(t))./sin(acos(t));

tD = zeros(S*3,2, Ng, M);
tN = zeros(S*3,2, Ng, M);

deltat = T/M;

uinc = zeros(2,Ng, M);
dnuinc = zeros(2,Ng, M);

dim_intf = 2*(2*Nmax +1);
lambdaD = zeros(S*3, dim_intf, M);
lambdaN = zeros(S*3, dim_intf, M);

for s = 1:S
gg = g(idx(s), :);

lambdaD((s-1)*3 + 1, :, :) = gg(1:dim_intf, :);
lambdaN((s-1)*3 + 1, :, :) = gg(dim_intf + 1:2*dim_intf, :);

lambdaD((s-1)*3 + 2, :, :) = gg(2*dim_intf + 1:3*dim_intf, :);
lambdaN((s-1)*3 + 2, :, :) = gg(3*dim_intf + 1:4*dim_intf, :);

lambdaD((s-1)*3 + 3, :, :) = gg(4*dim_intf + 1:5*dim_intf, :);
lambdaN((s-1)*3 + 3, :, :) = gg(5*dim_intf + 1:6*dim_intf, :);

clear gg;
end
%%
tic

for s = 1:S
for dom = 1:3


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



                           if length(J) == 1


                               tD((s-1)*3 + dom,jj, gg, m) = tD((s-1)*3 + dom,jj, gg, m) + (lambdaD((s-1)*3 + dom, counter + n + 1, m) * Un(n, xg(gg))./J);  

                               tN((s-1)*3 + dom,jj, gg, m) = tN((s-1)*3 + dom,jj, gg, m) + (lambdaN((s-1)*3 + dom, counter + n + 1, m) * Un(n, xg(gg))./J); 

                           else 

                               tD((s-1)*3 + dom,jj, gg, m) = tD((s-1)*3 + dom,jj, gg, m) + (lambdaD((s-1)*3 + dom, counter + n + 1, m) * Un(n, xg(gg))./J(gg));  

                               tN((s-1)*3 + dom,jj, gg, m) = tN((s-1)*3 + dom,jj, gg, m) + (lambdaN((s-1)*3 + dom, counter + n + 1, m) * Un(n, xg(gg))./J(gg)); 


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





%%
%%left half of the elipse;
x =@(t)a*cos(pi/2*(t+1)+pi-pi/2);
y = @(t) b*sin(pi/2*(t+1)+pi-pi/2);
dx = @(t) -a*pi/2*sin(pi/2*(t+1)+pi-pi/2);
dy =@(t) b*pi/2*cos(pi/2*(t+1)+pi-pi/2);
ddx = @(t) -a*pi*pi/4*cos(pi/2*(t+1)+pi-pi/2);
ddy = @(t) -b*pi*pi/4*sin(pi/2*(t+1)+pi-pi/2);

%%construction of the curve
c2 = Interface(x,y,dx,dy,ddx,ddy);

% c2 = {x, y, dx, dy, ddx, ddy, 1, 1, 28};

%%right half of the elipse;
x =@(t)a*cos(pi/2*(t+1)-pi/2);
y = @(t) b*sin(pi/2*(t+1)-pi/2);
dx = @(t) -a*pi/2*sin(pi/2*(t+1)-pi/2);
dy =@(t) b*pi/2*cos(pi/2*(t+1)-pi/2);
ddx = @(t) -a*pi*pi/4*cos(pi/2*(t+1)-pi/2);
ddy = @(t) -b*pi*pi/4*sin(pi/2*(t+1)-pi/2);

c1 = Interface(x,y,dx,dy,ddx,ddy);
% c1 = {x, y, dx, dy, ddx, ddy, 1, 1, 28};


%%separation between the two halves:
x=@(t) 0*ones(numel(t), 1);
y= @(t) -b*t;
dx =@(t) 0*ones(numel(t), 1);
dy =@(t) -b*1*ones(numel(t), 1);
ddx =@(t) 0*ones(numel(t), 1);
ddy =@(t) 0*ones(numel(t), 1);
%%we add a 0 at the end to indicate that this curve is not on the exterior
%domain
c3 =Interface(x,y,dx,dy,ddx,ddy,0);
% c3 = {x, y, dx, dy, ddx, ddy, 0, 1, 28};

%%separtion oriented in oposite direction.
c4 = c3.SetDirection(-1);
% c4 = {@(t)x(-t), @(t)y(-t), @(t) -dx(t), @(t) -dy(t), @(t)ddx(-t), @(t)ddy(-t), 0, 1, 28};


%%first domain Left domain.
%%NOT USE id =0!!!!

Interfz = cell(2,1);
Interfz{1} =c2;
Interfz{2} =c4;
k1=4;
d1 = Domain(1,0,-1.0*k1,Interfz);



%Second domain Right domain.
Interfz = cell(2,1);
Interfz{1} =c3;
Interfz{2} =c1;
k2 =8;
d2 =  Domain(2,0,-1.0*k2,Interfz);


% 
% 
% %upper half of the kite;
% x = @(t) a*(cos(pi/2*(t+1))+0.65*cos(2*(pi/2*(t+1)))-0.65);
% y = @(t) a*(1.5*sin(pi/2*(t+1)));
% dx = @(t) a*pi/2*(-sin(pi/2*(t+1))-1.3*sin(2*(pi/2*(t+1))));
% dy = @(t) a*pi/2*1.5*cos(pi/2*(t+1));
% ddx = @(t) a*pi/2*pi/2*(-cos(pi/2*(t+1))-2.6*cos(2*(pi/2*(t+1))));
% ddy = @(t) a*pi/2*pi/2*(-1.5*sin(pi/2*(t+1)));
% 
% %construction of the curve
% c2 = Interface(x,y,dx,dy,ddx,ddy);
% 
% %lower half of the kite;
% x = @(t) a*(cos(pi/2*(t+1)+pi)+0.65*cos(2*(pi/2*(t+1)+pi))-0.65);
% y = @(t) a*(1.5*sin(pi/2*(t+1)+pi));
% dx = @(t) a*pi/2*(-sin(pi/2*(t+1)+pi)-1.3*sin(2*(pi/2*(t+1)+pi)));
% dy = @(t) a*pi/2*1.5*cos(pi/2*(t+1)+pi);
% ddx = @(t) a*pi/2*pi/2*(-cos(pi/2*(t+1)+pi)-2.6*cos(2*(pi/2*(t+1)+pi)));
% ddy = @(t) a*pi/2*pi/2*(-1.5*sin(pi/2*(t+1)+pi));
% c1 = Interface(x,y,dx,dy,ddx,ddy);
% 
% %separation between the two halves:
% x=@(t) -a*t;
% y= @(t) 0*ones(numel(t),1);
% dx =@(t) -a*ones(numel(t),1);
% dy =@(t) 0*ones(numel(t),1);
% ddx =@(t) 0*ones(numel(t),1);
% ddy =@(t) 0*ones(numel(t),1);
% we add a 0 at the end to indicate that this curve is not on the exterior
% domain
% c3 =Interface(x,y,dx,dy,ddx,ddy,0);
% 
% %separtion oriented in oposite direction.
% c4 = c3.SetDirection(-1);
% 
% 
% %first domain Left domain.
% %NOT USE id =0!!!!
%  Interfz = cell(2,1);
%  Interfz{1} =c2;
%  Interfz{2} =c4;
%  k1=4;
%  d1 =  Domain(1,0,-1.0*k1,Interfz);
% 
% 
% 
%  %Second domain Right domain.
%  Interfz = cell(2,1);
%  Interfz{1} =c3;
%  Interfz{2} =c1;
%  k2 =8;
%  d2 =  Domain(2,0,-1.0*k2,Interfz);

%  % We save all the domains in one cell array
% Domains = cell(2,1);
% Domains{1} =d1;
% Domains{2} =d2;

%Exterior domain
c5 = c2.SetDirection(-1);
c6 = c1.SetDirection(-1);
Interfz{1} = c5;
Interfz{2} = c6;
d0 = Domain(0,0,-1.0*k2,Interfz);


[xg,wg] = lgwt(Ng,-1,1);
% xg = linspace(-1, 1, Ng+2);
% xg = xg(2:end-1).';
Un =@(n,t) sin((n+1)*acos(t))./sin(acos(t));


%% Evaluate the field

pause(2)
% 
Xlim = [-1 1];
Ylim = [-1 1];

xx = linspace(Xlim(1),Xlim(2),floor(abs(Xlim(2)-Xlim(1))*40));
yy = linspace(Ylim(1),Ylim(2),floor(abs(Ylim(2)-Ylim(1))*40));
[X,Y] = meshgrid(xx,yy);

Nx = numel(xx);
Ny = numel(yy);

pts = [reshape(X,Nx*Ny,1) reshape(Y,Nx*Ny,1)];


Npts = size(pts, 1);





% vert_1 = [x1(end:-1:1,:);[x1(1,1) h2];[x1(end,1) 100]];
[x1, y1] = d0.geo(xg, 1);
[x2, y2] = d0.geo(xg, 2);
dom0 = [x1 y1;x2 y2];

[x1, y1] = d1.geo(xg, 1);
[x2, y2] = d1.geo(xg, 2);
dom1 = [x1 y1;x2 y2];

[x1, y1] = d2.geo(xg, 1);
[x2, y2] = d2.geo(xg, 2);
dom2 = [x1 y1;x2 y2];



% in1 = inpolygon(pts(:,1),pts(:,2),vert_1(:,1),vert_1(:,2));
in1 = inpolygon(pts(:,1),pts(:,2),dom1(:,1),dom1(:,2));
in2 = inpolygon(pts(:,1),pts(:,2),dom2(:,1),dom2(:,2));
in0 = inpolygon(pts(:,1),pts(:,2),dom0(:,1),dom0(:,2));
in0 = ~in0;

%% POTENTIALS

disp('Potentials')

traces = cell(3*S,2);
lambda = R;
for s = 1:S
for dom = 1:3
dira = tD((s-1)*3 + dom, 1, :, :);
neua = tN((s-1)*3 + dom, 1, :, :);

dira = reshape(dira, Ng, M);
neua = reshape(neua, Ng, M);

dirb = tD((s-1)*3 + dom, 2, :, :);
neub = tN((s-1)*3 + dom, 2, :, :);

dirb = reshape(dirb, Ng, M);
neub = reshape(neub, Ng, M);

dir = zeros(2*Ng, M);
neu = zeros(2*Ng, M);

for n = 1:M
   dir(:, n) = [dira(:, n); dirb(:, n)];  
   neu(:, n) = [neua(:, n); neub(:, n)];
end

Lam = repmat(R.^(0:M-1),2*Ng,1);
dir = fft(Lam.*dir,[],2);
neu = fft(Lam.*neu,[],2);

traces{(s-1)*3 + dom, 1} = dir;
traces{(s-1)*3 + dom, 2} = neu;
end

end


ttraces = [];
for s = 1:S
    aux = [traces{(s-1)*3 + 1, 1};traces{(s-1)*3 + 1, 2}; ...
           traces{(s-1)*3 + 2, 1};traces{(s-1)*3 + 2, 2}; ...
           traces{(s-1)*3 + 3, 1};traces{(s-1)*3 + 3, 2}];
    ttraces = [ttraces; aux];
    clear aux;
end



%%
u_hlf = zeros(S*Npts,M);

idx=@(s) (s-1)*Npts+1:(s*Npts);
dim = 6*2*Ng;
idy=@(s) (s-1)*dim+1:(s*dim); 
hbar = parfor_progressbar(M,'Compute potentials...'); %create the progress bar 
tic

c = [c_0 c_1 c_2];
for n=1:M
    
    [P,Lambda]=eig(Am1-R*omega^(-(n-1))*B);
    Lambda=diag(Lambda)/dt;    
    
    gl=kron(inv(P),speye(dim))*ttraces(:,n);
      
    ul=zeros(S*Npts,1);
    
    
    for s = 1:S

    u_s = zeros(Npts, 1);
    
    
    k = 1i*Lambda(s)./c;


    mu = gl(idy(s));
    
    phi0 = mu(1:2*Ng); mu(1:2*Ng) = [];
    psi0 = mu(1:2*Ng); mu(1:2*Ng) = [];
    
    phi1 = mu(1:2*Ng); mu(1:2*Ng) = [];
    psi1 = mu(1:2*Ng); mu(1:2*Ng) = [];
    
    phi2 = mu(1:2*Ng); mu(1:2*Ng) = [];
    psi2 = mu(1:2*Ng); mu(1:2*Ng) = [];
    
    u0 = zeros(sum(in0), 1);    
    u1 = zeros(sum(in1), 1);
    u2 = zeros(sum(in2), 1);
    
    

    pot0 = compute_potentials(k(1),pts(in0,:),d0, Ng);

    pot1 = compute_potentials(k(2),pts(in1,:),d1, Ng);

    pot2 = compute_potentials(k(3),pts(in2,:),d2, Ng);
    
    

    u0 = pot0.SL*psi0-pot0.DL*phi0;

    u1 = pot1.SL*psi1-pot1.DL*phi1;
    
    u2 = pot2.SL*psi2-pot2.DL*phi2;
    
    
    % Putting everything together
    u_s(in0) = u0;
    u_s(in1) = u1;
    u_s(in2) = u2; 
    
    
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
    surf(X,Y,reshape(real(u(:, n)),Ny,Nx),'EdgeColor','black'); hold on;
    plot3(dom1(:,1),dom1(:,2),1*ones(size(dom1,1),1),'LineWidth',2,'color', 'k');
    plot3(dom2(:,1),dom2(:,2),1*ones(size(dom2,1),1),'LineWidth',2,'color', 'k');
    view(2);
    colormap(brewermap([],'*RdBu'))
    
%     colormap cool
%     colormap(brewermap([],'*RdGy'))
    caxis([-0.8 0.8])
    shading interp;
    hold off;
%     colorbar
%     axis equal
    axis off
%     axis([Xlim Ylim])
% %     
%     frm = getframe(gcf);
%     writeVideo(vid,frm)
%     

%     set(gca,'Xticklabel',[]) %to just get rid of the numbers but leave the ticks.
%     set(gca,'Yticklabel',[]) %to just get rid of the numbers but leave the ticks.

% drawnow;


% 

    
    count = count+1;
    pause(0.2)

    end

end


