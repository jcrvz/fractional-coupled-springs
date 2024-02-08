function [xt1, xt2] = Double_Mass_Spring_CF2(r,t,m1,k1,x1,v1,m2,k2,x2,v2)
%Double_Mass_Spring_CF fractional Solution of mechanical coupled 
%mass-spring system using the Caputo-Fabrizio fractional derivative
%  r gamma fractional order derivative
%  t time vector 
%  m1 mass 1
%  k1 constant of the spring 1
%  x1 initial position of m1, x1(t=0)
%  v1 initial velocity of the mass 1
%  m2 mass 2
%  k2 constant of the spring 2
%  x2 initial position of m2, x2(t=0)
%  v2 initial velocity of the mass 2
% created by Leonardo Martínez Jiménez 04/04/2020 

 %%% Z=sigma^(1-r) dimensionality factor (depends of the used value)
 Z=1.0;

 Br=(k2/m2+(k1+k2)/m1);
 Cr=4*k1*k2/(m1*m2);

 alfa2= Br/2-sqrt(Br^2-Cr)/2;    %%%alfa cuadrada
 beta2= Br/2+sqrt(Br^2-Cr)/2;    %%%beta cuadrada

 %%%%para x1
 a=(1-r)*(x1/m2+x2/m1)*k2*Z;
 b=(1-r)*(v1/m2+v2/m1)*k2*Z+r*(x1/m2+x2/m1)*k2*Z;
 c=r*(v1/m2+v2/m1)*k2*Z;
 %%%%parta x2
 a2=(1-r)*(k2*x1/m2+(k1+k2)*x2/m1)*Z;
 b2=(1-r)*(k2*v1/m2+(k1+k2)*v2/m1)*Z+r*(k2*x1/m2+(k1+k2)*x2/m1)*Z;
 c2=r*(k2*v1/m2+(k1+k2)*v2/m1)*Z;

 d=(1-r)*Z*alfa2;
 e=r*Z*alfa2;
 f=(1-r)*Z*beta2;
 g=r*Z*beta2;

 M1=[1 0 1 0;
     f 1 d 1;
     g f e d;
     0 g 0 e];
 ve1=[x1 v1+a b c];
 ve2=[x2 v2+a2 b2 c2];

 %%%%vector con los coeficientes 
 A1=inv(M1)*ve1';        %%% para x1  
 B1=inv(M1)*ve2';        %%% para x2 

 m=sqrt(r*Z*alfa2-(1-r)^2*Z*alfa2^2/4);
 n=sqrt(r*Z*beta2-(1-r)^2*Z*beta2^2/4);
 p=(1-r)*Z*alfa2/2;
 q=(1-r)*Z*beta2/2;

 %%%results
 xt1=A1(1)*exp(-p*t).*cos(m*t)+A1(2)/m*exp(-p*t).*sin(m*t)+A1(3)*exp(-q*t).*cos(n*t)+A1(4)/n*exp(-q*t).*sin(n*t);
 xt2=B1(1)*exp(-p*t).*cos(m*t)+B1(2)/m*exp(-p*t).*sin(m*t)+B1(3)*exp(-q*t).*cos(n*t)+B1(4)/n*exp(-q*t).*sin(n*t);
end