function solve_objtem
clc;clear;
h0=[0.1 0.5 0.1 0.1 0.1];
A = [1 1 1 1 1];
b = 1;
Aeq = [];beq = [];
lb = zeros(5,1);ub =ones(5,1);

options = optimset('Display','off','TolFun',1e-10);
options=optimset(options,'Algorithm','interior-point');
[x favel EXITFLAG]=fmincon(@obj,h0,A,b,Aeq,beq,lb,ub,@Non_cons,options)
end

%% 目标函数
function J=obj(h) 
rho1=1;rho2=1;A=100;N=6; 
y0=[10 0.5 10 11 16 ];
tspan=0:168;
for i=1:N;
    [tt,yb]=ode45(@state_NC,tspan,y0);  
    tspan=tspan+168;
    if i<N
          y0=[yb(end,1)*(1-h(i)) yb(end,2)*(1-h(i))^2 yb(end,3) yb(end,4) yb(end,5)*(1-h(i))];
    end 
end
J=A*rho1*sum(h)+rho2*yb(end,3);
end

%% 非线性约束条件
function [c ceq]=Non_cons(h)
delta=0.1;xi=150;N=6;
u=ODE(h);
b=length(u(:,1));
% g=Proj_Var(u,b);    
c=zeros(b,1);
for i=1:b-1;
    c(i)=u(i+1,1)+norminv(1-delta)*sqrt(u(i+1,2))-xi;
end

ceq=[];
end

% %% 投影
% function g=Proj_Var(u,b)
% x0=ones(b,2);
% A=[];b=[];Aeq=[];beq=[];VLB=[0 0];VUB=[];
% g=fmincon(@(x) Proj_err(x,u),x0,A,b,Aeq,beq,VLB,VUB,@VarConst);
% end
% 
% function v=Proj_err(x,u)
% v=sum(sum((x-u).^2,2));
% end
% 
% function [c,ceq]=VarConst(x)
% c=x(:,1).^2-x(:,2)+.001;   
% ceq=[];
% end
% 
%% 取出每段终端 N的一阶矩和二阶矩
function u=ODE(h)

u=[];
N=6;
y0=[10 0.5 10 11 16 ];
tspan=0:168;
N=6;

for i=1:N
      [tt yy]=ode45(@state_NC,tspan,y0);
      tspan=tspan+168;
      u(i,:)=[yy(end,1) yy(end,2)];
      if i<N
            y0=[yy(end,1)*(1-h(i)) yy(end,2)*(1-h(i))^2  yy(end,3) yy(end,4) yy(end,5)*(1-h(i))];           
      end
end
end


%%  状态方程
function dx =state_NC(t,x) 
% 1.EN;  2.EN^2;  3.EC;  4.EC^2;  5.ENC
alpha=0.0307;lambda=0.0115;eta=0.0000247;r=0.0026; 
dx=zeros(5,1);
  dx(1,1) = alpha + lambda*x(1) - eta*x(1)*x(3) - eta*x(5);
  dx(2,1) = alpha + 2*lambda*x(2) + lambda*x(1) + eta*x(1)*x(3) - 2*eta*x(1)*x(5) - 2*eta*x(2)*x(3) + eta*x(5);
  dx(3,1) = alpha + lambda*x(1) - r*x(3);
  dx(4,1) = alpha + 2*lambda*x(5) + lambda*x(1) - 2*r*x(4) + r*x(3);
  dx(5,1) = alpha + lambda*x(1) + lambda*x(5) + lambda*x(2) - eta*x(1)*x(4) - eta*x(3)*x(5) - r*x(5);
end