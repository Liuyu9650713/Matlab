
function g=Proj_Var(u,h)% transite from ODE solution to VarSpace with
% u 6 rows and 2 colms as input and return the same size of g��in VarSpace
% ���g�����㷽��������������,�������u���������������㡣
x0=ones(h,2);
A=[];b=[];Aeq=[];beq=[];VLB=[0 0];VUB=[];
g=fmincon(@(x) Proj_err(x,u),x0,A,b,Aeq,beq,VLB,VUB,@VarConst);
end
function v=Proj_err(x,u)% u is parameter or solution of ODE with six rows
v=sum(sum((x-u).^2,2));
end

function [c,ceq]=VarConst(x)% x,6 rows
c=x(:,1).^2-x(:,2)+.001;   
ceq=[];
end