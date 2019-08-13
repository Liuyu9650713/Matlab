function FigureMomentClosure
y0=[10 10 1 1 1];
tspan=[0 1000];
[t y]=ode45(@state_NC,tspan,y0,[]);
hold on
plot(t,y(:,2))
hold on
y(end,2)
end
% 1.EN;2.EC;3.EN^2;4.ENC;5.EC^2
% x(1)表示N的一阶矩，x(2)表示C的一阶矩，x(3)表示N的二阶矩，
% X(4)表示NC的混合矩，X(5)表示C的二阶矩，



function dx =state_NC(t,x) 
%alpha=0.03;lambda=0.012;eta=0.000025;r=0.003; 
alpha=0.0307;lambda=0.0115;eta=0.0000247;r=0.0026; 
dx=zeros(5,1);
  dx(1,1) = alpha - eta*x(4) + lambda*x(1);
  dx(2,1) = alpha + lambda*x(1) - r*x(2);
  dx(3,1) = alpha + eta*x(4) + 2*lambda*x(3) + x(1)*(2*alpha + lambda) - 2*eta*(2*x(1)*x(4) + x(2)*x(3) - 2*x(1)^2*x(2));
  dx(4,1) = alpha + x(1)*(alpha + lambda) + alpha*x(2) + lambda*x(3) - eta*(x(1)*x(5) + 2*x(2)*x(4) - 2*x(1)*x(2)^2) + x(4)*(lambda - r);
  dx(5,1) = alpha + lambda*x(1) + 2*lambda*x(4) - 2*r*x(5) + x(2)*(2*alpha + r);
end

% 1.EN;2.EC;3.EN^2;4.ENC;5.EC^2
% x(1)表示N的一阶矩，x(2)表示C的一阶矩，x(3)表示N的二阶矩，
% X(4)表示NC的混合矩，X(5)表示C的二阶矩，







% 
% 
% function dy =state_NC(t,y) %考虑只喷洒农药的矩方程
% alpha=0.0307;lambda=0.0115;eta=0.0000247;r=0.0026; 
% % y(1)表示N的一阶矩，y(2)表示N的二阶矩，y(3)表示C的一阶矩，y(4)表示C的二阶矩，y(5)表示NC的混合矩
% dy=zeros(5,1);
% dy=[alpha+lambda*y(1)-eta*y(1)*y(3)-eta*y(5);
%     alpha+2*lambda*y(2)+lambda*y(1)+eta*y(1)*y(3)-2*eta*y(1)*y(5)-2*eta*y(2)*y(3);
%     alpha+lambda*y(1)-r*y(3);
%     alpha+2*lambda*y(5)+lambda*y(1)-2*r*y(4)+r*y(3);
%     alpha+lambda*y(1)+lambda*y(5)+lambda*y(2)-eta*y(1)*y(4)-eta*y(3)*y(5)-eta*(y(1)*y(4)+2*y(3)*y(5)+y(1)*y(3)^2)-r*y(5)];
% end