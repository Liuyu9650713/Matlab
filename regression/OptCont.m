% a optimal proc by regression congstraint on parameter
function OptCont
A = [1 1 1 1 1  0 0 0 0 0; 0 0 0 0 0  1 1 1 1 1];
b = [1;1];
% the above means the sum of the proportions of epidemic is less than
% 1; and the sum of proportions of  steriles is less than 1 too.
% here theta=(mup1,mup2,...,mup6,mus1,mus2,...,mus6)' with 12 components.
% Atheta<=b
Aeq = [];beq = [];
%A=[];b=[];
% with no constaints of equality Aeq.theta=beq
lb = zeros(10,1);ub =ones(10,1);
% the lower bound and upper bound of theta is 0 and 1
ra1=rand(1,5);ra2=rand(1,5);theta1=ra1/sum(ra1)*rand;theta2=ra2/sum(ra2)*rand;
theta=[theta1 theta2];x0=theta;
%initial theta
options = optimset('Display','off','TolFun',1e-10);
options=optimset(options,'Algorithm','interior-point');
[x,fval,exitflag]=fmincon(@Obj,x0,A,b,Aeq,beq,lb,ub,@ContrF,options)
% where @ContrF is the nonlinear constraints c(x)<=0 and ceq(x)=0 which is
% defined below
%bar(x(1:5),0.2)

end

function d=Obj(x)
load Bgood
xx=[1;x(:)];
rhoc=1;
rhop=1;
rhos=4;
A=100;
Sbar=200;
d=sum(x(1:5))*rhop*A+sum(x(6:10))*rhos*Sbar+rhoc*exp(xx'*b6);
% the cost of pesticide and  steriles released and the last C left
end

function [c,ceq]=ContrF(x)
load B
%[mup(1)-mup(6),mus(1)-mus(6)] 
c(1)=exp([1 x(1) x(6)]*b1)-150;
c(2)=exp([1 x(1:2) x(6:7)]*b2)-150;
c(3)=exp([1 x(1:3) x(6:8)]*b3)-150;
c(4)=exp([1 x(1:4) x(6:9)]*b4)-150;
c(5)=exp([1 x(1:5) x(6:10)]*b5)-150;
ceq = [];
end