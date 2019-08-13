function [dx,magic]=mdyn4_fun(x)

%% created by closureDynamics(net,maxdeg=2,method='derMatch',funname='mdyn4_fun',symParameters='')

%% default parameters
  alpha = 0.0307;
  lambda = 0.0115;
  eta = 0.0000247;
  r = 0.0026;

if nargin==0
  % return variable names when nargin==0
  syms mu_N mu_C mu_N2 mu_N_C mu_C2
  dx={'mu_N';'mu_C';'mu_N2';'mu_N_C';'mu_C2'};
  magic=38033;  % random number used to check if correct function is being called
else
  % compute function
  dx(1,1) = alpha - eta*x(4) + lambda*x(1);
  dx(2,1) = alpha + lambda*x(1) - r*x(2);
  dx(3,1) = alpha + eta*x(4) + 2*lambda*x(3) + x(1)*(2*alpha + lambda) - (2*eta*x(3)*x(4)^2)/(x(1)^2*x(2));
  dx(4,1) = alpha + x(1)*(alpha + lambda) + alpha*x(2) + lambda*x(3) + x(4)*(lambda - r) - (eta*x(4)^2*x(5))/(x(1)*x(2)^2);
  dx(5,1) = alpha + lambda*x(1) + 2*lambda*x(4) - 2*r*x(5) + x(2)*(2*alpha + r);
end
