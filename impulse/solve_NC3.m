function solve_NC3
clc; clear;
global sigma;
d=3; beta=1.55;
sigma=10^2;

%T=20;
n=2;
N=2*n+2;
x0=zeros(1,2*n+2);
%x0(1:N-1)=[0.15 0.25 0.3 0.1 0.1];
x0(1:N-1)=[0.16 0.27 0.32 0.295 0.42];
x0(end)=0.1; %epsion

lb=[0.1*ones(1,N-1)   0]; %下界
ub=[5*ones(1,2*n+1)  5];

options=optimset('display','iter');
%options=optimset(options,'tolx',1e-8);
options=optimset(options,'GradObj','on');
AA=[1*ones(1,N-1) 0];
b=1;
[x,fval] = fmincon(@objfungrad,x0,AA,b,[],[],lb,ub,[],options)

end

function dy =state_NC1(t,y,tau)
alpha=0.03;lambda=0.012;eta=0.000025;r=0.003;rho1=1;rho2=1;A=100;delta=0.1;xi=150;

T=36;N=6;
tau=T/N;
dy=zeros(5,1);
dy=tau*[alpha+lambda*y(1)-eta*y(1)*y(3)-eta*y(5);
            alpha+2*lambda*y(2)+lambda*y(1)+eta*y(1)*y(3)-2*eta*y(1)*y(5)-2*eta*y(2)*y(3)+eta*y(5);
            alpha+lambda*y(1)-r*y(3);
            alpha+2*lambda*y(5)+lambda*y(1)-2*r*y(4)+r*y(3);
            alpha+lambda*y(1)+lambda*y(5)+lambda*y(2)-eta*y(1)*y(4)-eta*y(3)*y(5)-r*y(5)];
end
%function dy=costate_NC(t,y,tt,state1,state2,state3,state4,state5,[tau epsion])
function dy=costate_NC(t,y,tt,state1,state2,state3,state4,state5, epsion)
alpha=0.03;lambda=0.012;eta=0.000025;r=0.003;rho1=1;rho2=1;A=100;delta=0.1;xi=150;
d=3;
T=36;N=6;
tau=T/N;
dy=zeros(5,1);
y1=interp1(tt,state1,t);y2=interp1(tt,state2,t);y3=interp1(tt,state3,t);y4=interp1(tt,state4,t);y5=interp1(tt,state5,t);

dy=(-1)*(-1)*tau*[2*epsion^(-d)*max(y1+norminv(1-delta)*sqrt(y2-y1^2)-xi,0)*(1-norminv(1-delta)*y1/sqrt(y2-y1^2))+lambda*y(1)-eta*y(1)*y3+lambda*y(2)+eta*y(2)*y3-2*eta*y(2)*y5+lambda*y(3)+lambda*y(4)+lambda*y(5)-eta*y(5)*y4;
    2*epsion^(-d)*max(y1+norminv(1-delta)*sqrt(y2-y1^2)-xi,0)*norminv(1-delta)*1/(2*sqrt(y2-y1^2))+2*lambda*y(2)-2*eta*y(2)*y3+lambda*y(5);
    -eta*y(1)*y1+eta*y(2)*y1-2*eta*y(2)*y2-r*y(3)+r*y(4)-eta*y(5)*y5;
    -2*r*y(4)-eta*y(5)*y1;
    -eta*y(1)-2*eta*y(2)*y1+eta*y(2)+2*lambda*y(4)+lambda*y(5)-eta*y(5)*y3-r*y(5)];
end

function obj_J=obj(x)

global sigma;

alpha=0.03;lambda=0.012;eta=0.000025;r=0.003;rho1=1;rho2=1;A=100;delta=0.1;xi=150;
T=36;n=2;N=2*n+2;tau=T/N;

d=3; beta=1.55;

p=zeros(1,2*n+1);
p(1:2*n+1)=x(1:2*n+1);
epsion=x(end);
p(2*n+2)=0;

y0=[12 19 13 11 16];
tsw=zeros(1,2*n+2);
tsw(end)=T;

y=[];t=[];obj_conS=[];
tn=zeros(1,2*n+1);
tspan=[0:0.01:1];

for i=1:2*n+2;
    [tt,yy]=ode45(@state_NC1,tspan,y0,[],tau);
    tn(i)=length(tt);
    ya(i,:)=yy(1,:);yb(i,:)=yy(end,:);
    t=[t;tt];
    y=[y; yy];
    yyb=max(yy(:,1)+norminv(1-delta)*sqrt(yy(:,2)-yy(:,1).^2)-xi,0).^2;%列向量，约束违反函数中的约束 S(t)<xi1，
    obj_conSS=tau*trapz(tt,yyb);%数，约束违反函数中S(t)<xi1 的积分
    obj_conS=[obj_conS obj_conSS];%行向量
    y0=[yb(end,1)*(1-p(i)) yb(end,2)*(1-p(i))^2 yb(end,3)*(1-p(i)) yb(end,4) yb(end,5)];
    
end

con_viola=sum(obj_conS);%约束违反函数,数
obj_J=sum(A*rho1*p(1:2*n+1))+rho2*yb(end,3)+epsion^(-d)*con_viola+sigma*epsion^beta;

end

function grad_J=grad_NC(x)

global sigma;

alpha=0.03;lambda=0.012;eta=0.000025;r=0.003;rho1=1;rho2=1;A=100;delta=0.1;xi=150;
T=36;n=2;
N=2*n+2;tau=T/N;

d=3; beta=1.55;

p=zeros(1,2*n+1);
p(1:2*n+1)=x(1:2*n+1);
epsion=x(end);
p(2*n+2)=0;

tsw=zeros(1,2*n+3);
tsw(end)=T;

y0=[12 19 13 11 16];

y=[];t=[];obj_conS=[];
tn1=zeros(1,2*n+2);
tn2=zeros(1,2*n+1+1+1);
tspan=[0:0.1:1];

for i=1:2*n+2;
    [tt,yy]=ode45(@state_NC1,tspan,y0,[],tau);
    
    tn1(i)=length(tt);
    ya(i,:)=yy(1,:);yb(i,:)=yy(end,:);
    
    t=[t;tt];y=[y;yy];
    
    yyb=max(yy(:,1)+norminv(1-delta)*sqrt(yy(:,2)-yy(:,1).^2)-xi,0).^2;%列向量，约束违反函数中的约束 S(t)<xi1，
    obj_conSS=tau*trapz(tt,yyb);%数，约束违反函数中S(t)<xi1 的积分
    obj_conS=[obj_conS obj_conSS];%行向量
    
    if i<N
    y0=[yb(end,1)*(1-p(i))  yb(end,2)*(1-p(i))^2  yb(end,3)*(1-p(i))  yb(end,4)  yb(end,5)];
    end
end

tn2(2:end)=cumsum(tn1);

con_viola=sum(obj_conS);

grad_J=zeros(1,2*n+2);
grad_p=zeros(1,2*n+1);

lam_t=[];lamda_y=[];

tn3=zeros(1,2*n+1+1);
tn4=zeros(1,2*n+1+1+1);
lama=zeros(2*n+1,5);
lamb=zeros(2*n+1,5);

inv_tspan=[-1:0.1:0];

inv_lamdai0=[0 0 rho2 0 0];

for i=1:2*n+2;
    
    ti=t(tn2(N-i+1)+1:tn2(N-i+2));
    yi=y(tn2(N-i+1)+1:tn2(N-i+2),:);
    inv_ti=fliplr(-ti');inv_yi=fliplr(yi');
    [inv_tti,inv_lamdai]=ode45(@(t,y)costate_NC(t,y,inv_ti,inv_yi(1,:),inv_yi(2,:),inv_yi(3,:),inv_yi(4,:),inv_yi(5,:),epsion),inv_tspan,inv_lamdai0);
    tti=fliplr(-inv_tti');lamdai=fliplr(inv_lamdai');
    lam_t=[tti lam_t];lamda_y=[lamdai lamda_y];
    lama(N-i+1,:)=lamdai(:,1)'; lamb(N-i+1,:)=lamdai(:,end)';
    tn3(2*n+2-i+1)=length(tti);
    if i<N
        
        inv_lamdai0=[lama(N-i+1,1)*(1-p(N-i)) lama(N-i+1,2) lama(N-i+1,3) lama(N-i+1,4) lama(N-i+1,5)];
        
    end
end
tn4(2:end)=cumsum(tn3);

%下面把状态和协态进行标准化

for i=1:2*n+2;
    
    dt=0.02;  t_stan=0:dt:1;  %把协态和状态值统一在相同的离散点上
    
    %下面先把状态标准化
    ti=t(tn2(i)+1:tn2(i+1));%第i段上离散时间点，列向量
    yi=y(tn2(i)+1:tn2(i+1),:);%第i段上离散状态，第一列S的状态，第二列I的状态
    yi_stan=[interp1(ti,yi(:,1),t_stan) ;   interp1(ti,yi(:,2),t_stan); interp1(ti,yi(:,3),t_stan) ;  interp1(ti,yi(:,4),t_stan) ;  interp1(ti,yi(:,5),t_stan)];    %通过插值进行标准化，矩阵，第一行S的状态，第二行I的状态
    
    %把协态标准化
    tti=lam_t(tn4(i)+1:tn4(i+1));%第i段上的离散时间点 ，行向量
    lamdai=lamda_y(:,tn4(i)+1:tn4(i+1));%第i段上离散协态，第一行S的协态，第二行I的协态
    lamdai_stan=[interp1(tti,lamdai(1,:),t_stan);  interp1(tti,lamdai(2,:),t_stan);  interp1(tti,lamdai(3,:),t_stan);  interp1(tti,lamdai(4,:),t_stan);  interp1(tti,lamdai(5,:),t_stan)];  %通过插值进行标准化，矩阵，第一行S的状态，第二行I的状态
    
    if i<2*n+2;
        
        grad_p(i)=A*rho1-yb(i,1).*lama(i+1,1);
        
    end
    
    grad_J(end)=-d*epsion^(-d-1)*con_viola+sigma*beta*epsion^(beta-1);%对精度eps求导
    grad_J=[grad_p grad_J(end)];
    
end
end

function [f,G]=objfungrad(x)

f=obj(x);
G=grad_NC(x);

end









