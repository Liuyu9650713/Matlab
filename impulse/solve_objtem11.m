function solve_objtem11
h0=[0.1 0.1 0.1 0 0]; %Ͷ�ű����ĳ�ʼֵ
A = [1 1 1 1 1];
b = 1;
% sum of proportion must be less than 1
%A=[];b=[];
Aeq = [];beq = [];
lb = [zeros(5,1)];ub =[ones(5,1)];
% �����Ż�����������
options = optimset('Display','off','TolFun',1e-10);
options=optimset(options,'Algorithm','interior-point');
[x favel EXITFLAG]=fmincon(@obj,h0,A,b,Aeq,beq,lb,ub,@Non_cons,options)
end

function J=obj(h) % object function with Input h=[h1 h2 h3 h4 h5] Output cost of pesticides
rho1=1;rho2=1;A=100;N=6; y0=[10 119 10 120 59];
tspan=0:168;
for i=1:N-1;
    [tt,yb]=ode45(@state_NC,tspan,y0);  
    tspan=tspan+168;
    y0=[yb(end,1)*(1-h(i)) yb(end,2) yb(end,3)*(1-h(i))^2 yb(end,4)*(1-h(i)) yb(end,5)];
     %����ʹ��һʱ�̵ĳ�ʼ״̬�����ı䣨���壩
     % 1.EN;               2.EC;     3.EN^2;              4.ENC;             5.EC^2
     % x(1)��ʾN��һ�׾أ�x(2)��ʾC��һ�׾أ�x(3)��ʾN�Ķ��׾أ�
     % X(4)��ʾNC�Ļ�Ͼأ�X(5)��ʾC�Ķ��׾أ�
end
[tt,yb]=ode45(@state_NC,tspan,y0); 
J=A*rho1*sum(h)+rho2*yb(end,2); 
end

function [c ceq]=Non_cons(h)  % Nonlinear constraints with h=[h1 h2 h3 h4 h5]
delta=0.1;xi=150; u=ODE(h);  % u �ط�յĽ� 6 rows 2 cols
% g=Proj_Var(u);    % g ͶӰ�����ľ� 6 rows 2 cols
c=zeros(5,1);
for i=1:5
    c(i)=u(i+1,1)+norminv(1-delta)*sqrt(u(i+1,2)-u(i+1,1)^2)-xi;    % �������û�����ֳ�ѭ��
end
ceq=[];
end
% 
% function g=Proj_Var(u)% transite from ODE solution to VarSpace with
% % u 6 rows and 2 colms as input and return the same size of g��in VarSpace
% % ���g�����㷽��������������,�������u���������������㡣
% x0=ones(6,2);
% A=[];b=[];Aeq=[];beq=[];VLB=[0 0];VUB=[];
% g=fmincon(@(x) Proj_err(x,u),x0,A,b,Aeq,beq,VLB,VUB,@VarConst);
% end
% 
% function v=Proj_err(x,u)% u is parameter or solution of ODE with six rows
% v=sum(sum((x-u).^2,2));
% end
% 
% function [c,ceq]=VarConst(x)% x,6 rows
% c=x(:,1).^2-x(:,2)+.001;   
% ceq=[];
% end

function u=ODE(h)% h=[h1 h2 h3 h4 h5] is the control proportion of pesticide and
% u, a matirx with 6 rows 2 columns, is first two moment of N
% the first row is irrespective of what h is, and the other rows are
% determined by h
tspan=0:168;
N=6;
y0=[10 119 10 120 59];
[tt yy]=ode45(@state_NC,tspan,y0);
u(1,:)=[yy(end,1) yy(end,3)];%the first and the second monent of N
for i=2:N
      % 1.EN;               2.EC;     3.EN^2;                4.ENC;               5.EC^2
   y0=[yy(end,1)*(1-h(i-1)) yy(end,2) yy(end,3)*(1-h(i-1))^2 yy(end,4)*(1-h(i-1)) yy(end,5)];
   [tt yy]=ode45(@state_NC,tspan,y0);
   tspan=tspan+168;
   u(i,:)=[yy(end,1) yy(end,3)];
end
end

function dx =state_NC(t,x) 
%alpha=0.03;lambda=0.012;eta=0.000025;r=0.003; 
% 1.EN;2.EC;3.EN^2;4.ENC;5.EC^2
alpha=0.0307;lambda=0.0115;eta=0.0000247;r=0.0026; 
dx=zeros(5,1);
  dx(1,1) = alpha - eta*x(4) + lambda*x(1);
  dx(2,1) = alpha + lambda*x(1) - r*x(2);
  dx(3,1) = alpha + eta*x(4) + 2*lambda*x(3) + x(1)*(2*alpha + lambda) - 2*eta*(2*x(1)*x(4) + x(2)*x(3) - 2*x(1)^2*x(2));
  dx(4,1) = alpha + x(1)*(alpha + lambda) + alpha*x(2) + lambda*x(3) - eta*(x(1)*x(5) + 2*x(2)*x(4) - 2*x(1)*x(2)^2) + x(4)*(lambda - r);
  dx(5,1) = alpha + lambda*x(1) + 2*lambda*x(4) - 2*r*x(5) + x(2)*(2*alpha + r);
  
end