
function model()
clc;clear;
y0=[10 5 10 11 16 ];
t0=0;tfinal=168*6;
[T,Y]=ode45(@ODENC,[t0,tfinal],y0);

figure(1)
plot(T,Y(:,1),'R','linewidth',1.5)
xlabel('Time(hours)','FontSize',12);
ylabel('k10','FontSize',12);
set(gca,'linewidth',1.5);

figure(2)
plot(T,Y(:,2),'R','linewidth',1.5)
xlabel('Time(hours)','FontSize',12);
ylabel('k20','FontSize',12);
set(gca,'linewidth',1.5);

figure(3)
plot(T,Y(:,3),'R','linewidth',1.5)
xlabel('Time(hours)','FontSize',12);
ylabel('k01','FontSize',12);
set(gca,'linewidth',1.5);
% 
% figure(4)
% plot(T,Y(:,4),'R','linewidth',1.5)
% xlabel('Time(hours)','FontSize',12);
% ylabel('k02','FontSize',12);
% set(gca,'linewidth',1.5);
% 
% figure(5)
% plot(T,Y(:,5),'R','linewidth',1.5)
% xlabel('Time(hours)','FontSize',12);
% ylabel('k11','FontSize',12);
% set(gca,'linewidth',1.5);

end

function dydt=ODENC(t,y)
a=0.0307;lambda=0.0115;eta=0.247e-4;r=0.0026;rho1=1;rho2=1;m=100;de=0.1;xi=150;

dydt=[a+lambda*y(1)-eta*y(1)*y(3)-eta*y(5);
        a+2*lambda*y(2)+lambda*y(1)+eta*y(1)*y(3)-2*eta*y(1)*y(5)-2*eta*y(2)*y(3)+eta*y(5);
        a+lambda*y(1)-r*y(3);
        a+2*lambda*y(5)+lambda*y(1)-2*r*y(4)+r*y(3);
        a+lambda*y(1)+lambda*y(5)+lambda*y(2)-eta*y(1)*y(4)-eta*y(3)*y(5)-r*y(5)];
 end