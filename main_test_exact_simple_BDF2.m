
% clear 
close all
clc

%%
%-----------------------------------------
% set main parameters
%-----------------------------------------
eme = 1;
ene = 0;


% M = 10; % timesteps


emes =[25 50 100 200 400];
Mmax = 400;
% for Nmax = [20, eNes]


for M = [26]

disp('-----------------------------------------------')
% wavespeeds
Nmax = 20;
% M = emes(eme);

c_0 = 1.;
c_1 = 0.5;

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


%first domain Left domain.
%NOT USE id =0!!!!
 Interfz = cell(2,1);
 Interfz{1} =c2;
 Interfz{2} =c1;
 k1=4;
 d1 =  Domain(1,0,-1.0*k1,Interfz);

 % We save all the domains in one cell array
 Domains = cell(1,1);
 Domains{1} =d1;

%     M = 100;
T = 5; % 10
dt = T/M;
time = 0:dt:T;
% time = linspace(0, T, M);
tlag = 0.5;
% dirx = sqrt(0.5);
% diry = -sqrt(0.5);
dirx = 1;
diry = 0;
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


%%
disp(' ')
disp('TDMTF-solve');
tic
% hbar = parfor_progressbar(M/2+1,'Solving Linear Systems...'); %create the progress bar 
for l=0:M/2       %compute half the sequence

    if norm(solver.m_b(:,l+1)) > 1e-14
    k0 = p(R*omega^(-l))/(c_0 * T/M);
    k1 = p(R*omega^(-l))/(c_1 * T/M);

    %Set wave numbers for each subdomain
    solver.m_kext_real = real(k0);

    solver.m_kext_imag = imag(k0);

    solver.m_kext = k0;


    solver.m_DomainArray{1}.m_KwaveNumber = k0;
    solver.m_DomainArray{2}.m_KwaveNumber = k1;
    solver.ConstructDomain0();

    solver.setup(Nmax)
    rhs = solver.m_b(:,l+1);
    
    g(:,l+1)= solver.m_A \ rhs;  
    disp(l);

    end
    
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
Un =@(n,t) sin((n+1)*acos(t))./sin(acos(t));

tD = zeros(2,2, Ng, M+1);
tN = zeros(2,2, Ng, M+1);

deltat = T/M;

uinc = zeros(2,Ng, M+1);
dnuinc = zeros(2,Ng, M+1);

tD_ref = zeros(2,2,Ng, M+1);
tN_ref = zeros(2,2,Ng, M+1);

lambdaD = zeros(2, size(g, 1)/4, M+1);
lambdaN = zeros(2, size(g, 1)/4, M+1);


lambdaD(1, :, :) = g(1:size(g, 1)/4, :);
lambdaN(1, :, :) = g(size(g, 1)/4 + 1:size(g, 1)/2, :);

lambdaD(2, :, :) = g(size(g, 1)/2 + 1:3*size(g, 1)/4, :);
lambdaN(2, :, :) = g(3*size(g, 1)/4 + 1:end, :);
%%
tic
for dom = 1:2


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
                       onda = timedomain_test((m-1)*deltat, tlag, dirx, diry, c_0);
                       [dx, dy] = onda.evaluateGradient(z);

                       tN_ref(1,jj,:,m) = nx .* dx + ny .* dy;
                       tD_ref(1,jj,:, m) = onda.evaluate(z);
                       
                   elseif dom == 2

                       onda = timedomain_test((m-1)*deltat, tlag, dirx, diry, c_1);
                       [dx, dy] = onda.evaluateGradient(z);

                       tN_ref(2,jj,:,m) = nx .* dx + ny .* dy;
                       tD_ref(2,jj,:, m) = onda.evaluate(z);
                             
                   end


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



%left half of the elipse;
x =@(t)a*cos(pi/2*(t+1)+pi-pi/2);
y = @(t) b*sin(pi/2*(t+1)+pi-pi/2);
dx = @(t) -a*pi/2*sin(pi/2*(t+1)+pi-pi/2);
dy =@(t) b*pi/2*cos(pi/2*(t+1)+pi-pi/2);
ddx = @(t) -a*pi*pi/4*cos(pi/2*(t+1)+pi-pi/2);
ddy = @(t) -b*pi*pi/4*sin(pi/2*(t+1)+pi-pi/2);

%construction of the curve
c2 = Interface(x,y,dx,dy,ddx,ddy);

% c2 = {x, y, dx, dy, ddx, ddy, 1, 1, 28};

%right half of the elipse;
x =@(t)a*cos(pi/2*(t+1)-pi/2);
y = @(t) b*sin(pi/2*(t+1)-pi/2);
dx = @(t) -a*pi/2*sin(pi/2*(t+1)-pi/2);
dy =@(t) b*pi/2*cos(pi/2*(t+1)-pi/2);
ddx = @(t) -a*pi*pi/4*cos(pi/2*(t+1)-pi/2);
ddy = @(t) -b*pi*pi/4*sin(pi/2*(t+1)-pi/2);

c1 = Interface(x,y,dx,dy,ddx,ddy);
% c1 = {x, y, dx, dy, ddx, ddy, 1, 1, 28};



%first domain Left domain.
%NOT USE id =0!!!!

Interfz = cell(2,1);
Interfz{1} =c2;
Interfz{2} =c1;
k1=4;
d1 = Domain(1,0,-1.0*k1,Interfz);



%Exterior domain
c5 = c2.SetDirection(-1);
c6 = c1.SetDirection(-1);
Interfz{1} = c5;
Interfz{2} = c6;
d0 = Domain(0,0,-1.0*k1,Interfz);


k_hlf = 1i*p(R*omega.^(-(0:floor(M/2))))/dt;
[xg,wg] = lgwt(Ng,-1,1);
Un =@(n,t) sin((n+1)*acos(t))./sin(acos(t));


%% Evaluate the field

pause(2)
% 
Xlim = [-2 2];
Ylim = [-2 2];

xx = linspace(Xlim(1),Xlim(2),floor(abs(Xlim(2)-Xlim(1))*2.5));
yy = linspace(Ylim(1),Ylim(2),floor(abs(Ylim(2)-Ylim(1))*2.5));
[X,Y] = meshgrid(xx,yy);

Nx = numel(xx);
Ny = numel(yy);

% pts = [reshape(X,Nx*Ny,1) reshape(Y,Nx*Ny,1)];
% pts = [2 0;-0.5 0;0.5 0];
pts = [2 0;0 0];


% vert_1 = [x1(end:-1:1,:);[x1(1,1) h2];[x1(end,1) 100]];
[x1, y1] = d0.geo(xg, 1);
[x2, y2] = d0.geo(xg, 2);
dom0 = [x1 y1;x2 y2];

[x1, y1] = d1.geo(xg, 1);
[x2, y2] = d1.geo(xg, 2);
dom1 = [x1 y1;x2 y2];



% in1 = inpolygon(pts(:,1),pts(:,2),vert_1(:,1),vert_1(:,2));
in1 = inpolygon(pts(:,1),pts(:,2),dom1(:,1),dom1(:,2));
in0 = inpolygon(pts(:,1),pts(:,2),dom0(:,1),dom0(:,2));
in0 = ~in0;
%%
% Incident field



disp('Potentials')

traces = cell(2,2);
lambda = R;
for dom = 1:2
dira = tD(dom, 1, :, :);
neua = tN(dom, 1, :, :);

dira = reshape(dira, Ng, M+1);
neua = reshape(neua, Ng, M+1);

dirb = tD(dom, 2, :, :);
neub = tN(dom, 2, :, :);

dirb = reshape(dirb, Ng, M+1);
neub = reshape(neub, Ng, M+1);

dir = zeros(2*Ng, M+1);
neu = zeros(2*Ng, M+1);

for n = 1:M+1
   dir(:, n) = [dira(:, n); dirb(:, n)];  
   neu(:, n) = [neua(:, n); neub(:, n)];
end

Lam = repmat(lambda.^(0:M),2*Ng,1);
dir = fft(Lam.*dir,[],2);
neu = fft(Lam.*neu,[],2);

traces{dom, 1} = dir;
traces{dom, 2} = neu;
end


clear lambdaD lambdaN


u_hlf = zeros(Nx*Ny,M/2+1);

tic

c = [c_0 c_1];
for n=1:M/2+1

    k = k_hlf(n)./c;


    
    phi0 = traces{1, 1}(:, n);
    psi0 = traces{1, 2}(:, n);
    
    phi1 = traces{2, 1}(:, n);
    psi1 = traces{2, 2}(:, n);
    
    
    u0 = zeros(sum(in0), 1);    
    u1 = zeros(sum(in1), 1);
    
    

    pot0 = compute_potentials(k(1),pts(in0,:),d0, Ng);

    pot1 = compute_potentials(k(2),pts(in1,:),d1, Ng);
    
    

    u0 = pot0.SL*psi0-pot0.DL*phi0;

    u1 = pot1.SL*psi1-pot1.DL*phi1;
    % Putting everything together
    u_hlf(in0, n) = u0;
    u_hlf(in1, n) = u1;
toc
end
toc

u_hlf(:,M+2-(1:floor(M/2)))=conj(u_hlf(:,2:floor(M/2)+1));  % mirror the hermitian seq
u=real(ifft(u_hlf,[],2));                            % w
u=bsxfun(@times,u,R.^(-(0:M)));                % g


u0 = u(in0, :);

u1 = u(in1, :);

clear u dir neu traces



if eme == 0
    u0_ref = u0;
    u1_ref = u1;
    tD_ref = tD;
    tN_ref = tN;
    
else
disp('Errors')



errorD(eme, 1) = traceErrorDom(Ai, 0*tD_ref, tD,0, 1)*deltat;
errorD(eme, 2) = traceErrorDom(Ai, tD_ref, tD,0, 2)
% errorD(eme, 4) = traceError(Ai, tD_ref, tD,0, 2, 2)
% errorD(eme, 5) = traceError(Ai, tD_ref, tD,0, 3, 1)
% errorD(eme, 6) = traceError(Ai, tD_ref, tD,0, 3, 2)




errorN(eme, 1) = traceErrorDom(Ai, 0*tN_ref, tN,1, 1)*deltat;
errorN(eme, 2) = traceErrorDom(Ai, tN_ref, -tN,1, 2)
% errorN(eme, 4) = traceErrorDom(Ai, tN_ref, -tN,1, 2, 2)
% errorN(eme, 5) = traceErrorDom(Ai, tN_ref, -tN,1, 3, 1)
% errorN(eme, 6) = traceErrorDom(Ai, tN_ref, -tN,1, 3, 2)

% 
% z = pts(in0, 1) + 1i*pts(in0, 2);
% u0_ref = zeros(1, M+1);
% for m = 1:M+1
% onda = timedomain_test((m-1)*deltat, tlag, dirx, diry, c_0);
% u0_ref(m) = onda.evaluate(z);
% end


z = pts(in1, 1) + 1i*pts(in1, 2);
u1_ref = zeros(1, M+1);
for m = 1:M+1
onda = timedomain_test((m-1)*deltat, tlag, dirx, diry, c_1);
u1_ref(m) = onda.evaluate(z);
end


% ind = 1:Mmax/M:Mmax+1;
errorU(eme, 1) = norm(u0)*deltat; 
errorU(eme, 2) = norm(u1_ref - u1)/norm(u1_ref)



save('test-simple.mat', 'errorD', 'errorN', 'errorU');

end
clearvars inc inc2;



save(['exact-solution-simple-',num2str(M),'.mat']);


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


% 
% save('test.mat', 'errorD', 'errorN', 'errorU');
% % mail_results('Resultados TDMTF', date);
% %%
% % emes = [50 100 200 400];
% % % % subplot(2, 1, 1);
loglog(emes, errorD(:, 1), '-s');hold on;
loglog(emes, errorD(:, 2), '-^');
loglog(emes, 11*emes.^(-3/2), '--k');
loglog(emes, 0.75e4*emes.^(-2), '-.k');
legend({'dom0', 'dom1','Order 3/2', 'Order 2'},...
        'Location', 'southwest');hold off;
% 

% 
% saveas(gcf, 'errorD.png');
% close all
% % % % % % subplot(2, 1, 2);    
loglog(emes, errorN(:, 1), '-s');hold on;
loglog(emes, errorN(:, 2), '-^');
loglog(emes, 55*emes.^(-1), '--k');
loglog(emes, 0.75e2*emes.^(-3/2), '-.k');
legend({'dom0', 'dom1','Order 1', 'Order 3/2'},...
        'Location', 'southwest');
% 
% saveas(gcf, 'errorN.png');
% close all
% % 
% 
% % 
loglog(emes, errorU(:, 1), '-s');hold on;
loglog(emes, errorU(:, 2), '-^');
loglog(emes, 11*emes.^(-3/2), '--k');
loglog(emes, 0.75e4*emes.^(-2), '-.k');
legend({'0', '1','Order 3/2', 'Order 2'},...
        'Location', 'southwest');hold off;

