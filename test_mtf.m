clear 
close all
clc

%%
%-----------------------------------------
% set main parameters
%-----------------------------------------
eme = 0;


enes = [80 2 3 5 7 9 11 13 15 25 40];

freqs = [1-1i 2 4 8];

for freq = 1
    
%     
% s0 = freqs(freq);
% s1 = 5*2*s0;
% s2 = 5*4*s0;
% 

% s0 = freqs(freq);
% s1 = 5*2*s0;
% s2 = 25*4*s0;

s0 = freqs(freq);
s1 = 5*2*s0;
s2 = 25*4*s0;

ene = 0;

for Nmax = enes

disp('-----------------------------------------------')
disp(Nmax)



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




kwave = 1i*s0;
theta = 0;
inc{1} = plane_wave(theta, kwave);



kwave_real = 0;
kwave_imag = -4*pi;


%%
solver = solverMtf28(1,kwave_real,kwave_imag,inc,Domains);


solver.setupb_freq(Nmax)
disp('Lado derecho')

        
        



%%
tic
disp(' ')
disp('MTF-solve');

k0 = s0;
k1 = s1;
k2 = s2;

%Set wave numbers for each subdomain
solver.m_kext_real = real(k0);

solver.m_kext_imag = imag(k0);

solver.m_kext = k0;


solver.m_DomainArray{1}.m_KwaveNumber = k0;
solver.m_DomainArray{2}.m_KwaveNumber = k1;
solver.m_DomainArray{3}.m_KwaveNumber = k2;
solver.ConstructDomain0();

solver.setup(Nmax)


rhs = solver.m_b(:,1);

g = solver.m_A \ rhs;  

toc

disp(' ')





%%

% disp('Computing Traces');
% Ng = 2*max(solver.m_N{1});

Ng = 101;
% [xg,wg] = lgwt(Ng,-1,1);
xg = linspace(-1, 1, Ng+2);
xg = xg(2:end-1).';
Un =@(n,t) sin((n+1)*acos(t))./sin(acos(t));

tD = zeros(3,2, Ng);
tN = zeros(3,2, Ng);



uinc = zeros(2,Ng);
dnuinc = zeros(2,Ng);


lambdaD = zeros(3, size(g, 1)/6);
lambdaN = zeros(3, size(g, 1)/6);


lambdaD(1, :) = g(1:size(g, 1)/6);
lambdaN(1, :) = g(size(g, 1)/6 + 1:size(g, 1)/3);

lambdaD(2, :) = g(size(g, 1)/3 + 1:size(g, 1)/2);
lambdaN(2, :) = g(size(g, 1)/2 + 1:2*size(g, 1)/3);

lambdaD(3, :) = g(2*size(g, 1)/3 + 1:5*size(g, 1)/6);
lambdaN(3, :) = g(5*size(g, 1)/6 + 1:size(g, 1));

for dom = 1:3
    counter = 0;
     NIntrfz0 = solver.m_DomainArray{dom}.m_NumInfefaces;
     for jj=1:NIntrfz0
       J = solver.m_DomainArray{dom}.J(xg, jj);
       for n = 0:2*solver.m_N{dom}(jj)
        for gg=1:Ng
           if length(J) == 1
               tD(dom,jj, gg) = tD(dom,jj, gg) + (lambdaD(dom, counter + n + 1) * Un(n, xg(gg))./J);  
               tN(dom,jj, gg) = tN(dom,jj, gg) + (lambdaN(dom, counter + n + 1) * Un(n, xg(gg))./J); 
           else 
               tD(dom,jj, gg) = tD(dom,jj, gg) + (lambdaD(dom, counter + n + 1) * Un(n, xg(gg))./J(gg));  
               tN(dom,jj, gg) = tN(dom,jj, gg) + (lambdaN(dom, counter + n + 1) * Un(n, xg(gg))./J(gg)); 
           end
        end
       end
     end
       counter = counter + 2*solver.m_N{dom}(jj) + 1;
end


%% Errors



load('Ai.mat')


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


%separation between the two halves:
x=@(t) 0*ones(numel(t), 1);
y= @(t) -b*t;
dx =@(t) 0*ones(numel(t), 1);
dy =@(t) -b*1*ones(numel(t), 1);
ddx =@(t) 0*ones(numel(t), 1);
ddy =@(t) 0*ones(numel(t), 1);
% we add a 0 at the end to indicate that this curve is not on the exterior
% domain
c3 =Interface(x,y,dx,dy,ddx,ddy,0);
% c3 = {x, y, dx, dy, ddx, ddy, 0, 1, 28};

%separtion oriented in oposite direction.
c4 = c3.SetDirection(-1);
% c4 = {@(t)x(-t), @(t)y(-t), @(t) -dx(t), @(t) -dy(t), @(t)ddx(-t), @(t)ddy(-t), 0, 1, 28};


%first domain Left domain.
%NOT USE id =0!!!!

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

%Exterior domain
c5 = c2.SetDirection(-1);
c6 = c1.SetDirection(-1);
Interfz{1} = c5;
Interfz{2} = c6;
d0 = Domain(0,0,-1.0*k2,Interfz);


% [xg,wg] = lgwt(Ng,-1,1);
xg = linspace(-1, 1, Ng+2);
xg = xg(2:end-1).';
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
pts = [1 0;-0.25 0;0.25 0];

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
%%
% Incident field



% disp('Potentials')

traces = cell(3,2);

for dom = 1:3
dira = tD(dom, 1, :);
neua = tN(dom, 1, :);

dira = reshape(dira, Ng, 1);
neua = reshape(neua, Ng, 1);

dirb = tD(dom, 2, :);
neub = tN(dom, 2, :);

dirb = reshape(dirb, Ng, 1);
neub = reshape(neub, Ng, 1);


dir = [dira; dirb];  
neu = [neua; neub];

traces{dom, 1} = dir;
traces{dom, 2} = neu;
end



u_hlf = zeros(Nx*Ny);


k = 1i*[s0 s1 s2];



phi0 = traces{1, 1};
psi0 = traces{1, 2};

phi1 = traces{2, 1};
psi1 = traces{2, 2};

phi2 = traces{3, 1};
psi2 = traces{3, 2};

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
u_hlf(in0) = u0;
u_hlf(in1) = u1;
u_hlf(in2) = u2; 





u0 = u_hlf(in0);

u1 = u_hlf(in1);

u2 = u_hlf(in2);

clear u dir neu traces



disp('Errors')

if ene == 0
    gref = g;
    tD_ref = tD;
    tN_ref = tN;
    u0_ref = u0;

    u1_ref = u1;

    u2_ref = u2;
    
    lambdaN_ref = lambdaN;
    lambdaD_ref = lambdaD;
    
else
errorD(ene, 1) = traceErrorDom_freq(Ai, gref, g,0, 1,1);
errorD(ene, 2) = traceErrorDom_freq(Ai, gref, g,0, 1,2);
errorD(ene, 3) = traceErrorDom_freq(Ai, gref, g,0, 2,1);
errorD(ene, 4) = traceErrorDom_freq(Ai, gref, g,0, 2,2);
errorD(ene, 5) = traceErrorDom_freq(Ai, gref, g,0, 3,1);
errorD(ene, 6) = traceErrorDom_freq(Ai, gref, g,0, 3,2);


errorN(ene, 1) = traceErrorDom_freq(Ai, gref, g,1, 1,1);
errorN(ene, 2) = traceErrorDom_freq(Ai, gref, g,1, 1,2);
errorN(ene, 3) = traceErrorDom_freq(Ai, gref,g,1, 2,1);
errorN(ene, 4) = traceErrorDom_freq(Ai, gref, g,1, 2,2);
errorN(ene, 5) = traceErrorDom_freq(Ai, gref, g,1, 3,1);
errorN(ene, 6) = traceErrorDom_freq(Ai, gref,g,1, 3,2);

errorU(ene, 1) = norm(u0_ref - u0)/norm(u0_ref); 
errorU(ene, 2) = norm(u1_ref - u1)/norm(u1_ref); 
errorU(ene, 3) = norm(u2_ref - u2)/norm(u2_ref); 


% save(['test_freq-',num2str(freq), '.mat'], 'errorD', 'errorN', 'errorU');
% mail_results(['MTF - ', num2str(Nmax)], compose([num2str(Nmax), '\n', ...
%                         'errorD: ', num2str(errorD(end, :)),'\n', ...
%                         'errorN: ', num2str(errorN(end, :)),'\n', ...
%                         'errorU: ', num2str(errorU(end, :))]), 'test_freq.mat');

clearvars inc inc2;


end
ene = ene +1;
% save(['solution_freq-',num2str(Nmax),'.mat']);


end


% mail_results('MTF', num2str(freq), ['test_freq-',num2str(freq), '.mat']);

end


% semilogy(2*[2 3 5 7 9 11 13 15]+1, errorD(:, :), '-o');hold on;
% semilogy(2*[2 3 5 7 9 11 13 15]+1, errorN(:, :), '-x');
% semilogy(2*[2 3 5 7 9 11 13 15]+1, errorU(:, :), '-s');


% loglog(2*[2 3 5 7 9 11 13 15]+1, errorD(:, :), '-o');hold on;
% loglog(2*[2 3 5 7 9 11 13 15]+1, errorN(:, :), '-x');
% loglog(2*[2 3 5 7 9 11 13 15]+1, errorU(:, :), '-s');
loglog(2*[2 3 5 7 9 11 13 15 25 40]+1, errorD(:, :), '-ro');hold on;
loglog(2*[2 3 5 7 9 11 13 15 25 40]+1, errorN(:, :), '-bx');
loglog(2*[2 3 5 7 9 11 13 15 25 40]+1, errorU(:, :), '-s');

emes = 2*[2 3 5 7 9 11 13 15 25 40]+1;

slopesN = zeros(numel(emes)-1, 6);

for nn = 1:(numel(emes)-1)
   slopesN(nn, 1) = -log(errorN(nn+1, 1)/errorN(nn,1))/log(emes(nn+1)/emes(nn));  
   slopesN(nn, 2) = -log(errorN(nn+1, 2)/errorN(nn,2))/log(emes(nn+1)/emes(nn));  
   slopesN(nn, 3) = -log(errorN(nn+1, 3)/errorN(nn,3))/log(emes(nn+1)/emes(nn));  
   slopesN(nn, 4) = -log(errorN(nn+1, 4)/errorN(nn,4))/log(emes(nn+1)/emes(nn));  
   slopesN(nn, 5) = -log(errorN(nn+1, 5)/errorN(nn,5))/log(emes(nn+1)/emes(nn));  
   slopesN(nn, 6) = -log(errorN(nn+1, 6)/errorN(nn,6))/log(emes(nn+1)/emes(nn));  
end



slopesD = zeros(numel(emes)-1, 6);

for nn = 1:(numel(emes)-1)
   slopesD(nn, 1) = -log(errorD(nn+1, 1)/errorD(nn,1))/log(emes(nn+1)/emes(nn));  
   slopesD(nn, 2) = -log(errorD(nn+1, 2)/errorD(nn,2))/log(emes(nn+1)/emes(nn));  
   slopesD(nn, 3) = -log(errorD(nn+1, 3)/errorD(nn,3))/log(emes(nn+1)/emes(nn));  
   slopesD(nn, 4) = -log(errorD(nn+1, 4)/errorD(nn,4))/log(emes(nn+1)/emes(nn));  
   slopesD(nn, 5) = -log(errorD(nn+1, 5)/errorD(nn,5))/log(emes(nn+1)/emes(nn));  
   slopesD(nn, 6) = -log(errorD(nn+1, 6)/errorD(nn,6))/log(emes(nn+1)/emes(nn));  
end


slopesU = zeros(numel(emes)-1, 3);

for nn = 1:(numel(emes)-1)
   slopesU(nn, 1) = -log(errorU(nn+1, 1)/errorU(nn,1))/log(emes(nn+1)/emes(nn));  
   slopesU(nn, 2) = -log(errorU(nn+1, 2)/errorU(nn,2))/log(emes(nn+1)/emes(nn));  
   slopesU(nn, 3) = -log(errorU(nn+1, 3)/errorU(nn,3))/log(emes(nn+1)/emes(nn));   
end



save('test-freq4.mat', 'errorD', 'errorN', 'errorU', 'emes')




