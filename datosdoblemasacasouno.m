%%% Doble Masa (Dr. Rosales)
clear
set(0,'DefaultAxesFontSize',14);
set(0,'defaultLineLineWidth',1.5);
csvread caso1.csv;
t = ans(1:595,1);
l1 = ans(1:599,2);
l2 = ans(1:599,3);
x1 = l1(4:598);
x2 = l2(4:598);
m1=0.12723;
m2=0.12781;
k1=10.57372949;
k2=9.792093254;
%ofd=.95;atenuacion rapida
ofd=0.99675;%sin atenuacion
%ofd=0.9977;%si funciona
al = [ofd,ofd,ofd,ofd];
tp = 5.39106e-44;
A = -(tp^(1-ofd))*(k1+k2)/m1; B = (tp^(1-ofd))* k2/m1;
C = (tp^(1-ofd))*k2/m2; D = -(tp^(1-ofd))*k2/m2;
param = [A,B,C,D];
lam = [0 , 1, 0, 1];
f_fun=@(r,y,par)[...
    y(2);...
    A*y(1)+B*y(3);...
    y(4);...
    C*y(1)+D*y(3)];

J_fun=@(r,y,par)[...
    0,1,0,0;...
    A,0,B,0;...
    0,0,0,1;...
    C,0,D,0];

r0=0;R=20;
y0=[0.0175;0;0.0175;0];
h = 0.0001;
%[r, y] = MT_FDE_PI1_Ex(al,lam,f_fun,r0,R,y0,h,param);
[r, y] = mt_fde_pi1_im(al,lam,f_fun,J_fun,r0,R,y0,h,param);
z = y(1,:);
zv = y(2,:);
w = y(3,:);
wv = y(4,:);

%%%%  CAPUTO-FABRIZZIO
%%%%condiciones iniciales
tcf=0:0.01:20;
x1cf=0.035; x2cf=0.035;            %%%%posiciones iniciales
v1=0; v2=0;           %%%%velocidades iniciales 
ordfrac=0.9977;%%%%orden derivación fraccional
[xt1, xt2] = Double_Mass_Spring_CF2(ordfrac,tcf,m1,k1,x1cf,v1,m2,k2,x2cf,v2);

%%%% Graficas cartesianas x1 vs x2
% %%%%%
% subplot(1,2,1)
% plot(x1,x2,'b',z,w,'r')
% title('$x_1(t)\ vs\ x_2(t)\ \textrm{by Caputo}$','Interpreter','latex')
% xlabel('$x_{1}(t)\ [mts]$','Interpreter','latex')
% ylabel('$x_{2}(t)\ [mts]$','Interpreter','latex')
% legend({'\textrm{Data}\ $(x_1,x_2)$','$(x_1(t),x_2(t))\ \textrm{Caputo}$'},'Interpreter','latex','Location','northwest','FontSize',14)
% subplot(1,2,2)
% plot(x1,x2,'b',xt1,xt2,'r')
% title('$x_1(t)\ vs\ x_2(t)\ \textrm{by Caputo-Fabrizio}$','Interpreter','latex')
% xlabel('$x_{1}(t)\ [mts]$','Interpreter','latex')
% ylabel('$x_{2}(t)\ [mts]$','Interpreter','latex')
% legend({'\textrm{Data}\ $(x_1,x_2)$','$(x_1(t),x_2(t))\ \textrm{Caputo-Fabrizio}$'},'Interpreter','latex','Location','northwest','FontSize',14)
%%%%%%%%%%%%%
%%%%%%%%%%%%%%



subplot(2,1,1)
plot(t,x1,'r',r,z,'b',tcf,xt1,'g','LineWidth',2)
title('$\textrm{Case 1.}\ m_1=0.12723\ kg,\ m_2=0.12781\ kg,\ k_1=10.57372949\ N/m,k_2=9.792093254\ N/m$','Interpreter','latex')
xlabel('$t\ [seg]$','Interpreter','latex')
ylabel('$x_{1}(t)\ [mts]$','Interpreter','latex')
legend({'\textrm{Data}\ $x_1$','$x_1(t)\ \textrm{Caputo}$','$x_1(t)\ \textrm{C Fabrizzio}$'},'Interpreter','latex','Location','northeast','FontSize',14)
axis([0 20 -0.04 0.04])
subplot(2,1,2)
%plot(t,x2,'r')
plot(t,x2,'r',r,w,'b',tcf,xt2,'g','LineWidth',2)
legend({'\textrm{Data}\ $x_2$','$x_2(t)\ \textrm{Caputo}$','$x_1(t)\ \textrm{C Fabrizzio}$'},'Interpreter','latex','Location','northeast','FontSize',14)
xlabel('$t\ [seg]$','Interpreter','latex')
ylabel('$x_{2}(t)\ [mts]$','Interpreter','latex')
axis([0 20 -0.065 0.065])

%%%%%%%% Distancia discreta de Frechet
% for k=1:100;
% taux(k)=t(6*k-5);
% x1aux(k)=x1(6*k-5);
% x2aux(k)=x2(6*k-5);
% end
% for k=1:100;
% raux(k)=r(2000*k+1);
% zaux(k)=z(2000*k+1);
% waux(k)=w(2000*k+1);
% end
% %%%% Distancia de los datos a Caputo x1
% P1=[taux' x1aux'];
% Q1=[raux' zaux'];
% [cm1, cSq1] = DiscreteFrechetDist(P1,Q1);
% 
% % plot result
% figure(2)
% plot(Q1(:,1),Q1(:,2),'o-r','linewidth',3,'markerfacecolor','r')
% hold on
% plot(P1(:,1),P1(:,2),'o-b','linewidth',3,'markerfacecolor','b')
% title(['Discrete Frechet Distance of curves P and Q: ' num2str(cm1)])
% legend({'$x_1(t)\ \textrm{by Caputo}$','$\textrm{Data}\ x_1(t)$'},'Interpreter','latex','location','best')
% xlabel('$t\ [seg]$','Interpreter','latex')
% ylabel('$x_{1}(t)\ [mts]$','Interpreter','latex')
% line([2 cm1+2],[0.035 0.035],'color','m','linewidth',2)
% text(2+cm1+0.1,0.035,'dFD length')
% for i=1:length(cSq1)
%   line([P1(cSq1(i,1),1) Q1(cSq1(i,2),1)],...
%        [P1(cSq1(i,1),2) Q1(cSq1(i,2),2)],...
%        'color',[0 0 0]+(i/length(cSq1)/1.35));
% end
% 
% %%%% Distancia de los datos a Caputo x2
% P2=[taux' x2aux'];
% Q2=[raux' waux'];
% [cm2, cSq2] = DiscreteFrechetDist(P2,Q2);
% 
% % plot result
% figure(3)
% plot(Q2(:,1),Q2(:,2),'o-r','linewidth',3,'markerfacecolor','r')
% hold on
% plot(P2(:,1),P2(:,2),'o-b','linewidth',3,'markerfacecolor','b')
% title(['Discrete Frechet Distance of curves P and Q: ' num2str(cm2)])
% legend({'$x_2(t)\ \textrm{by Caputo}$','$\textrm{Data}\ x_2(t)$'},'Interpreter','latex','location','best')
% xlabel('$t\ [seg]$','Interpreter','latex')
% ylabel('$x_{2}(t)\ [mts]$','Interpreter','latex')
% line([2 cm2+2],[0.05 0.05],'color','m','linewidth',2)
% text(2+cm2+0.1,0.05,'dFD length')
% for i=1:length(cSq2)
%   line([P2(cSq2(i,1),1) Q2(cSq2(i,2),1)],...
%        [P2(cSq2(i,1),2) Q2(cSq2(i,2),2)],...
%        'color',[0 0 0]+(i/length(cSq2)/1.35));
% end
% 
% 
% for k=1:100;
% tcfaux(k)=tcf(20*k+1);
% xt1aux(k)=xt1(20*k+1);
% xt2aux(k)=xt2(20*k+1);
% end
% 
% P1=[taux' x1aux'];
% Q3=[tcfaux' xt1aux'];
% [cm3, cSq3] = DiscreteFrechetDist(P1,Q3);
% 
% % plot result
% figure(4)
% plot(Q3(:,1),Q3(:,2),'o-r','linewidth',3,'markerfacecolor','r')
% hold on
% plot(P1(:,1),P1(:,2),'o-b','linewidth',3,'markerfacecolor','b')
% title(['Discrete Frechet Distance of curves P and Q: ' num2str(cm3)])
% legend({'$x_1(t)\ \textrm{by Caputo-Fabrizio}$','$\textrm{Data}\ x_1(t)$'},'Interpreter','latex','location','best')
% xlabel('$t\ [seg]$','Interpreter','latex')
% ylabel('$x_{1}(t)\ [mts]$','Interpreter','latex')
% line([2 cm3+2],[0.035 0.035],'color','m','linewidth',2)
% text(2+cm3+0.1,0.035,'dFD length')
% for i=1:length(cSq3)
%   line([P1(cSq3(i,1),1) Q3(cSq3(i,2),1)],...
%        [P1(cSq3(i,1),2) Q3(cSq3(i,2),2)],...
%        'color',[0 0 0]+(i/length(cSq3)/1.35));
% end
% 
% P2=[taux' x2aux'];
% Q4=[tcfaux' xt2aux'];
% [cm4, cSq4] = DiscreteFrechetDist(P2,Q4);
% 
% % plot result
% figure(5)
% plot(Q4(:,1),Q4(:,2),'o-r','linewidth',3,'markerfacecolor','r')
% hold on
% plot(P2(:,1),P2(:,2),'o-b','linewidth',3,'markerfacecolor','b')
% title(['Discrete Frechet Distance of curves P and Q: ' num2str(cm4)])
% legend({'$x_2(t)\ \textrm{by Caputo-Fabrizio}$','$\textrm{Data}\ x_2(t)$'},'Interpreter','latex','location','best')
% xlabel('$t\ [seg]$','Interpreter','latex')
% ylabel('$x_{2}(t)\ [mts]$','Interpreter','latex')
% line([2 cm4+2],[0.05 0.05],'color','m','linewidth',2)
% text(2+cm4+0.1,0.05,'dFD length')
% for i=1:length(cSq4)
%   line([P2(cSq4(i,1),1) Q4(cSq4(i,2),1)],...
%        [P2(cSq4(i,1),2) Q4(cSq4(i,2),2)],...
%        'color',[0 0 0]+(i/length(cSq4)/1.35));
% end