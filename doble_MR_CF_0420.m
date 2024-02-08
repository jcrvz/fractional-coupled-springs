clear all
close all
clc

set(0,'DefaultAxesFontSize',20);
set(0,'defaultLineLineWidth',1.5);

%%%base de tiempo
t=0:0.01:60;

m1=1; m2=1;            %%%%masas del sistema
k1=6; k2=4;            %%%%constantes del resorte

%%%%condiciones iniciales
x1=-1; x2=2;            %%%%posiciones iniciales
v1=0; v2=0;           %%%%velocidades iniciales 

%%%parametros fraccionarios
r=0.970;  %%%% orden fraccionario

[xt1, xt2] = Double_Mass_Spring_CF2(r,t,m1,k1,x1,v1,m2,k2,x2,v2);
figure, plot(t,xt1,t,xt2,'r');
xlabel('$t$  $[s]$','Interpreter','latex' );
ylabel('$x$ $[m]$','Interpreter','latex');
legend({'$x_1(t)$','$x_2(t)$'},'Interpreter','latex');
%axis([0 60 -0.75 0.75])


