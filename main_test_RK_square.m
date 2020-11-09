%% RUNGE KUTTA CONVOLUTION QUADRATURE

clear 
% load 'square2-solution-RK1-100.mat'
close all
clc


%% TEST

eme = 0;
ene = 0;
Mmax = 400;
emes =[50 100 200];

%% Runge-Kutta parameters
Ark=[5/12 -1/12; 3/4 1/4];   % RADAU IIa
% Ark=[1/6 -1/3 1/6; 1/6 5/12 -1/12; 1/6 2/3 1/6]; % LOBATTO IIIc
bb=Ark(end,:);
S=size(Ark,1); RK = S;
Am1 = inv(Ark);
crk=Ark*ones(S,1);
B = (Am1*ones(size(Ark,1),1))*[zeros(1,S-1),1];
Nmax = 10;
for M = [Mmax emes]
%% PDE parameters
disp('-----------------------------------------------')
% wavespeeds

% M = emes(eme);

c_0 = 1.;
c_1 = 0.75;
c_2 = 0.5;
c_3 = 0.25;
c_4 = 0.125;



%% DOMAINs

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

%%
 
T = 10;
dt = T/M;
time = 0:dt:T;
% time = linspace(0, T, M);
tlag = 1;
dirx = sqrt(0.5);
diry = -sqrt(0.5);

omega = exp(2*pi*1i/(M));
R = eps^(0.5/(M));

%% Right hand side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inc{1} = td_test(T-dt, M-1, 1,tlag,dirx, diry, c_0, c_1, crk(1)*dt, 1);
inc2{1} = td_test(T-dt, M-1, 1,tlag,dirx, diry, c_0, c_1, crk(2)*dt, 1);
% inc3{1} = td_test(T-dt, M-1, 1,tlag,dirx, diry, c_0, c_1, crk(3)*dt, 1);

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

% if M == 200
%     load 'sol-ram.mat'
%     integers = 21:M-1;
% else
%     integers = 0:M-1;
% end
integers = 0:M-1;
%%
disp(' ')
disp('TDMTF-solve');
tic
% hbar = parfor_progressbar(M,'Solving Linear Systems...'); %create the progress bar 
for l=integers

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
    
%     hbar.iterate(1);
%     save('sol-ram.mat','g');
    toc
%     pause(20)
end
toc
% close(hbar);
disp(['Number of timesteps: ', num2str(M)]);

disp(' ')

g=real(ifft(g,[],2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%

g=bsxfun(@times,g,R.^(-(0:M-1)));
disp('Computing Traces');



Ng = 101;
[xg,wg] = lgwt(Ng,-1,1);
% xg = linspace(-1, 1, Ng+2);
% xg = xg(2:end-1).';
Un =@(n,t) sin((n+1)*acos(t))./sin(acos(t));



tD0 = zeros(1*S,8, Ng, M);
tN0 = zeros(1*S,8, Ng, M);


tD = zeros(4*S,4, Ng, M);
tN = zeros(4*S,4, Ng, M);

deltat = T/M;

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
                   
%                    if dom == 1 && s == S
%                        onda = timedomain_test(m*deltat, tlag, dirx, diry, c_0);
%                        [dx, dy] = onda.evaluateGradient(z);
% 
%                        tN_ref(1,jj,:,m) = nx .* dx + ny .* dy;
%                        tD_ref(1,jj,:, m) = onda.evaluate(z);
%                        
%                    elseif dom == 2
%                        onda = timedomain_test(m*deltat, tlag, dirx, diry, c_1);
%                        [dx, dy] = onda.evaluateGradient(z);
% 
%                        tN_ref(2,jj,:,m) = nx .* dx + ny .* dy;
%                        tD_ref(2,jj,:, m) = onda.evaluate(z);
%                    elseif dom == 3
%                        onda = timedomain_test(m*deltat, tlag, dirx, diry, c_2);
%                        [dx, dy] = onda.evaluateGradient(z);
% 
%                        tN_ref(3,jj,:,m) = nx .* dx + ny .* dy;
%                        tD_ref(3,jj,:, m) = onda.evaluate(z);
%                    end
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

% pts = [reshape(X,Nx*Ny,1) reshape(Y,Nx*Ny,1)];

% Evaluation points
% pts = [1 0;-0.25 0;0.25 0];

pts = [b/2 a/2;-b/2 a/2;-b/2 -a/2;b/2 -a/2; 2*b 2*a];
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
% hbar = parfor_progressbar(M,'Compute potentials...'); %create the progress bar 
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
    
%     hbar.iterate(1);
toc
end
toc
% close(hbar);

Lam = repmat(R.^(0:M-1),S*Npts,1);

u = Lam.^(-1).*ifft(u_hlf,[],2);

u = u((S-1)*Npts+1:end,:);



u0 = [zeros(sum(in0), 1) u(in0, :)];
u1 = [zeros(sum(in1), 1) u(in1, :)];
u2 = [zeros(sum(in2), 1) u(in2, :)];
u3 = [zeros(sum(in3), 1) u(in3, :)];
u4 = [zeros(sum(in4), 1) u(in4, :)];

clear u dir neu traces ttraces





%% ERRORS


load 'Ai.mat'

disp('Errors')
tD0 = tD0(S, :, :, :);
tN0 = tN0(S, :, :, :);
% 
tD = tD(5:8, :, :, :);
tN = tN(5:8, :, :, :);
% 
% tD = tD(9:12, :, :, :);
% tN = tN(9:12, :, :, :);

auxD = zeros(4,4, Ng, M+1);
auxN = zeros(4,4, Ng, M+1);

auxD(:,:,:, 2:M+1) = tD;
auxN(:,:,:, 2:M+1) = tN;

tD = auxD;
tN = auxN;


auxD0 = zeros(1,8, Ng, M+1);
auxN0 = zeros(1,8, Ng, M+1);

auxD0(:,:,:, 2:M+1) = tD0;
auxN0(:,:,:, 2:M+1) = tN0;

tD0 = auxD0;
tN0 = auxN0;

clear auxD auxN auxD0 auxN0;

if eme == 0
    u0_ref = u0;
    u1_ref = u1;
    u2_ref = u2;
    u3_ref = u3;
    u4_ref = u4;
    tD_ref = tD;
    tN_ref = tN;
    tD0_ref = tD0;
    tN0_ref = tN0;
    
    
else
disp('Errors')

load 'Ai_square.mat'
cte = 101*8;

A0 = Ai(1:2*cte, 1:2*cte); A0 = A0(1:size(A0,1)/2, size(A0, 1)/2+1:end);
A1 = Ai(2*cte + (1:cte),2*cte + (1:cte)); A1 = A1(1:size(A1,1)/2, size(A1, 1)/2+1:end);
A2 = Ai(3*cte + (1:cte),3*cte + (1:cte)); A2 = A2(1:size(A2,1)/2, size(A2, 1)/2+1:end);
A3 = Ai(4*cte + (1:cte),4*cte + (1:cte)); A3 = A3(1:size(A3,1)/2, size(A3, 1)/2+1:end);
A4 = Ai(5*cte + (1:cte),5*cte + (1:cte)); A4 = A4(1:size(A4,1)/2, size(A4, 1)/2+1:end);

errorD(eme, 1) = traceErrorDom_square(A0, tD0_ref, tD0,0, 1);
errorD(eme, 2) = traceErrorDom_square(A1, tD_ref, tD,0, 2);
errorD(eme, 3) = traceErrorDom_square(A2, tD_ref, tD,0, 3);
errorD(eme, 4) = traceErrorDom_square(A3, tD_ref, tD,0, 4);
errorD(eme, 5) = traceErrorDom_square(A4, tD_ref, tD,0, 5)


errorN(eme, 1) = traceErrorDom_square(A0, tN0_ref, tN0,1, 1);
errorN(eme, 2) = traceErrorDom_square(A1, tN_ref, tN,1, 2);
errorN(eme, 3) = traceErrorDom_square(A2, tN_ref, tN,1, 3);
errorN(eme, 4) = traceErrorDom_square(A3, tN_ref, tN,1, 4);
errorN(eme, 5) = traceErrorDom_square(A4, tN_ref, tN,1, 5)

ind = 1:Mmax/M:Mmax+1;
errorU(eme, 1) = norm(u0_ref(ind) - u0)/norm(u0_ref(ind)); 
errorU(eme, 2) = norm(u1_ref(ind) - u1)/norm(u1_ref(ind)); 
errorU(eme, 3) = norm(u2_ref(ind) - u2)/norm(u2_ref(ind)); 
errorU(eme, 4) = norm(u3_ref(ind) - u3)/norm(u3_ref(ind));
errorU(eme, 5) = norm(u4_ref(ind) - u4)/norm(u4_ref(ind))


% save('testRK.mat', 'errorD', 'errorN', 'errorU');
% mail_results(['TDMTF - Radau ', num2str(M)], [num2str(M), '                          ', ...
%                         'errorD: ', num2str(errorD(end, :)),'                          ', ...
%                         'errorN: ', num2str(errorN(end, :)),'                          ', ...
%                         'errorU: ', num2str(errorU(end, :))], 'testRK.mat');

end

% 
% 
% errorD(eme, 1) = traceErrorDom(Ai, 0*tD_ref, tD,0, 1)*deltat
% errorD(eme, 2) = traceErrorDom(Ai, tD_ref, tD,0, 2)
% errorD(eme, 3) = traceErrorDom(Ai, tD_ref, tD,0, 3)
% 
% errorN(eme, 1) = traceErrorDom(Ai, 0*tN_ref, -tN,1, 1)*deltat
% errorN(eme, 2) = traceErrorDom(Ai, tN_ref, -tN,1, 2)
% errorN(eme, 3) = traceErrorDom(Ai, tN_ref, -tN,1, 3)
% 
% % 
% % z = pts(in0, 1) + 1i*pts(in0, 2);
% % u0_ref = zeros(1, M+1);
% % for m = 1:M
% % onda = timedomain_test(m*deltat, tlag, dirx, diry, c_0);
% % u0_ref(m) = onda.evaluate(z);
% % end
% 
% 
% z = pts(in1, 1) + 1i*pts(in1, 2);
% u1_ref = zeros(1, M);
% for m = 1:M
% onda = timedomain_test(m*deltat, tlag, dirx, diry, c_1);
% u1_ref(m) = onda.evaluate(z);
% end
% 
% 
% 
% z = pts(in2, 1) + 1i*pts(in2, 2);
% u2_ref = zeros(1, M);
% for m = 1:M
% onda = timedomain_test(m*deltat, tlag, dirx, diry, c_1);
% u2_ref(m) = onda.evaluate(z);
% end
% 
% errorU(eme, 1) = norm(u0)*deltat
% errorU(eme, 2) = norm(u1_ref - u1)/norm(u1_ref) 
% errorU(eme, 3) = norm(u2_ref - u2)/norm(u2_ref) 
% 
% 
% 
% save('testRK.mat', 'errorD', 'errorN', 'errorU');
% mail_results(['TDMTF - Runge-Kutta - ', num2str(M)], [num2str(M), '\n', ...
%                         'errorD: ', num2str(errorD(end, :)),'\n', ...
%                         'errorN: ', num2str(errorN(end, :)),'\n', ...
%                         'errorU: ', num2str(errorU(end, :))], 'testRK.mat');

clearvars inc inc2 inc3;



save(['square2-solution-RK1-',num2str(M),'.mat']);

eme = eme +1;


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



save('test_RK1_square2.mat', 'ErrorD', 'ErrorN', 'ErrorU');
% mail_results('Resultados TDMTF', date);
%%
% emes = [50 100 200 400];
% % % subplot(2, 1, 1);
loglog(emes, errorD(:, 1), '-s');hold on;
loglog(emes, errorD(:, 2), '-^');
loglog(emes, errorD(:, 3), '-d');
loglog(emes, errorD(:, 4), '-o');
loglog(emes, errorD(:, 5), '-v');
loglog(emes, 11e2*emes.^(-3), '--k');
loglog(emes, 0.75e5*emes.^(-4), '-.k');
legend({'dom0', 'dom1', 'dom2','dom3','dom4','Order 3', 'Order 2'},...
        'Location', 'southwest');hold off;
% 


saveas(gcf, 'errorD-RK.png');
close all
% % % % % subplot(2, 1, 2);    
loglog(emes, errorN(:, 1), '-s');hold on;
loglog(emes, errorN(:, 2), '-^');
loglog(emes, errorN(:, 3), '-d');
loglog(emes, errorN(:, 4), '-o');
loglog(emes, errorN(:, 5), '-v');
loglog(emes, 55*emes.^(-1), '--k');
loglog(emes, 0.75e2*emes.^(-3/2), '-.k');
legend({'dom0', 'dom1', 'dom2','Order 1/2', 'Order 1'},...
        'Location', 'southwest');

saveas(gcf, 'errorN-RK.png');
close all
% 
% 
% % 
loglog(emes, errorU(:, 1), '-s');hold on;
loglog(emes, errorU(:, 2), '-^');
loglog(emes, errorU(:, 3), '-d');
loglog(emes, errorU(:, 4), '-o');
loglog(emes, errorU(:, 5), '-v');
loglog(emes, 11*emes.^(-3/2), '--k');
loglog(emes, 0.75e4*emes.^(-3), '-.k');
legend({'0', '1', '2','Order 1', 'Order 2'},... 
'Location', 'southwest');hold off;


saveas(gcf, 'errorU-RK.png');
close all

% mail_results('TDMTF-RungeKutta', 'Results', {'errorD-RK.png', 'errorN-RK.png', 'errorU-RK.png'})



