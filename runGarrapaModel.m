% Settings for plot
set(0,'DefaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1.5);

% Temporal array
t = linspace(0, 10, 100);

% Model parameters
m1 = 0.12723;       % kg
m2 = 0.12781;       % kg
k1 = 10.57372949;   % N/m
k2 = 9.792093254;   % N/m

% Fractional order
ofd = 0.99675;

% 
al = ofd * ones(1, 4);

% Model's parameters
A = -(k1+k2)/m1; 
B = k2/m1;
C = k2/m2; 
D = -k2/m2;

param = [A, B, C, D];
lam = [0 , 1, 0, 1];

f_fun = @(r, y, par) ...
    [y(2); 
     A*y(1)+B*y(3);
     y(4); 
     C*y(1)+D*y(3)];

J_fun = @(r,y,par) ...
    [0, 1, 0, 0;
     A, 0, B, 0;
     0, 0, 0, 1;
     C, 0, D, 0];

r0 = 0;
R = 20;
y0 = [0.0175; 0; 0.0175; 0];
h = 0.0001;

[r, y] = mt_fde_pi1_im(al, lam,f_fun,J_fun,r0,R,y0,h,param);
z = y(1,:);
zv = y(2,:);
w = y(3,:);
wv = y(4,:);

subplot(2,1,1)
plot(r,z,'b','LineWidth',2)
title('$\textrm{Case 1.}\ m_1=0.12723\ kg,\ m_2=0.12781\ kg,\ k_1=10.57372949\ N/m,k_2=9.792093254\ N/m$','Interpreter','latex')
xlabel('$t\ [seg]$','Interpreter','latex')
ylabel('$x_{1}(t)\ [mts]$','Interpreter','latex')
legend({'$x_1(t)\ \textrm{Caputo}$'},'Interpreter','latex','Location','northeast','FontSize',14)
axis([0 20 -0.04 0.04])