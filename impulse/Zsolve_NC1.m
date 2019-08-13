function Zsolve_NC1
clc ;clear;
N=5;
x0=zeros(1,N);
x0(1:N)=[0.1 0.2 0.25 0.3 0.2];


lb=[0*ones(1,N) 0];
ub=[1*ones(1,N) 0.1];

options = optimset('display','iter');%迭代的每个值显示出来
options = optimset(options,'GradObj','on');  %用梯度迭代

A=ones(1,N);
b=1;

[x,fval]=fmincon(@Z_objgrade1,x0,A,b,[],[],lb,ub,@Z_Nostrain,options)
end

%% 状态
function dy =Z_state1(t,y,tau)
a=0.0307;lambda=0.0115;eta=0.247e-4;r=0.0026;rho1=1;rho2=1;m=100;de=0.1;xi=150;

dy=zeros(5,1);
dy=tau*[a+lambda*y(1)-eta*y(1)*y(3)-eta*y(5);
        a+2*lambda*y(2)+lambda*y(1)+eta*y(1)*y(3)-2*eta*y(1)*y(5)-2*eta*y(2)*y(3)+eta*y(5);
        a+lambda*y(1)-r*y(3);
        a+2*lambda*y(5)+lambda*y(1)-2*r*y(4)+r*y(3);
        a+lambda*y(1)+lambda*y(5)+lambda*y(2)-eta*y(1)*y(4)-eta*y(3)*y(5)-r*y(5)];
end

%% 协态方程
function dy=Z_costate1(t,y,tt,state1,state2,state3,state4,state5)
a=0.0307;lambda=0.0115;eta=0.247e-4;r=0.0026;rho1=1;rho2=1;m=100;

T=36;N=5;
tau=T/(N+1);

dy=zeros(5,1);
y1=interp1(tt,state1,t);
y2=interp1(tt,state2,t);
y3=interp1(tt,state3,t);
y4=interp1(tt,state4,t);
y5=interp1(tt,state5,t);

dy=(-1)*tau*[lambda*y(1)-eta*y(1)*y3+lambda*y(2)+eta*y(2)*y3-2*eta*y(2)*y5+lambda*y(3)+lambda*y(4)+lambda*y(5)-eta*y(5)*y4;
             2*lambda*y(2)-2*eta*y(2)*y3+lambda*y(5);
             -eta*y(1)*y1+eta*y(2)*y1-2*eta*y(2)*y2-r*y(3)+r*y(4)-eta*y(5)*y5;
             -2*r*y(4)-eta*y(5)*y1;
             -eta*y(1)-2*eta*y(2)*y1+eta*y(2)+2*lambda*y(4)+lambda*y(5)-eta*y(5)*y3-r*y(5)];
end

%% 目标函数和梯度
function [obj_f,grad_G]=Z_objgrade1(x)
a=0.0307;lambda=0.0115;eta=0.247e-4;r=0.0026;rho1=1;rho2=1;m=100;de=0.1;xi=150;s=100;

T=36;N=5;tau=T/(N+1);

p=zeros(1,N);
p(1:N)=x(1:N);


y0=[4 19 3 11 16 ];

% 求目标
ya=zeros(N+1,5);yb=zeros(N+1,5);%把每一段初始值和结束时刻的状态值记录下来 ya:初始  yb:结束
y=[];t=[];%空矩阵开始默认是零矩阵，对应求解出来的时间点，状态
tn1=zeros(1,N+1); %每一段上离散点的个数取长之后构成的行向量，如tn1=[3,5,3,...]第一段上3个离散点，第2段上5个离散点，第3段上3个离散点...
tn2=zeros(1,N+2);%2*n+2个时间点上对应的离散点数构成的行向量,如tn2=[0,3,8,11,...]从第2个点开始等于tn1的累加求和
obj_conS=[]; 

 for i=1:N+1;%第一部分：调用状态函数，求解每一点的状态
         tspan=[i-1:0.1:i];
         [tt,yy] = ode45(@Z_state1,tspan,y0,[],tau);%tt是第i段上离散时间点的列向量，yy状态矩阵，第一列是S的数量，第二列是I的数量，形式yy(?,:)，其中？表示第几行，对应离散时刻，：表示第几个状态
         tn1(i)=length(tt);%每一段变换到[0 1]上，第i段离散时间点取长即点的数目
         
         ya(i,:)=yy(1,:); yb(i,:)=yy(end,:); %（行向量）把状态在0和1时刻的值记录下来,':'代表第几个状态,'i'代表第i个时间段，1代表0时刻，end代表1时刻，最后累计变成矩阵，第一列S，第二列I
         t=[t;tt]; %每一次把离散时间点的列向量累计起来，构成列向量t
         y=[y; yy]; %  y is a length(tt)*2 vector; 把每次的状态累计起来,第一列S，第二列I
         
         
         if i<N+1
             y0=[yb(i,1)*(1-p(i)) yb(i,2)*(1-p(i))^2 yb(i,3) yb(i,4) yb(i,5)*(1-p(i))];               
         end
 end

 tn2(2:end)=cumsum(tn1); %第2个时间点t1到最后一个时间点tN等于每一段的离散点的个数累加求和,如tn2=[0,3,8,11,...，20,23,25]从第2个点开始等于tn1的累加求和

obj_f=sum(m*rho1*p(1:N))+rho2*yb(end,3);

%求梯度
grad_G=zeros(1,N+1);%定义梯度规则
grad_p=zeros(1,N);
grad_epsion=zeros(1,1);
  
 %下面开始调用协态函数求解协态 
lam_t=[]; lamda_y=[];%空矩阵开始默认是零矩阵，对应求解出来的时间点，协态
lama=zeros(N+1,5);
lamb=zeros(N+1,5);
tn3=zeros(1,N+1); %每一段上离散点的个数取长之后构成的行向量，同tn1
tn4=zeros(1,N+2);%2*n+2个时间点上对应的离散点数构成的行向量,同tn2

inv_lamdai0=[0 0 rho2 0 0];%协态的在最后一段时间点1上的初始值

 for i=1:N+1;
     
     %i对应第N-i+1个时间段先求第5段，然后第4段，...
     
     ti=t(tn2(N+1-i+1)+1:tn2(N+1-i+2));%每一段上的时间点构成的列向量,（对应t里的元素），如果最后一段是两个时间点0.2,0.9，则t1=[0.2 0.9]',...,t5=[0.3 0.4 0.8]',倒求
     yi=y(tn2(N+1-i+1)+1:tn2(N+1-i+2),:);%把每一段上每个点的状态记录下来 （对应y里的元素），‘：’代表状态1,状态2,第一列代表S，第2列代表I，倒求
     inv_tspan=[N+1-i+1:-0.1:N+1-i];

     [inv_tti,inv_lamdai]=ode45(@(t,y)Z_costate1(t,y,ti,yi(:,1),yi(:,2),yi(:,3),yi(:,4),yi(:,5)),inv_tspan, inv_lamdai0);
     %inv_tti是第N-i+1段上离散时间点列向量，inv_lamdai是协态的矩阵，第一列是S的协态，第二列是I的协态，形式inv_lamdai(?,:)，其中？表示第几个时刻，：表示第几个协态

     tti=fliplr(inv_tti'); lamdai=fliplr(inv_lamdai'); %时间，协态翻转
     lama(N+1-i+1,:)=lamdai(:,1)'; lamb(N+1-i+1,:)=lamdai(:,end)';
     tn3(N+1-i+1)=length(tti);%每一段上对应协态离散化的时间点取长,tn3倒取
     
     lam_t=[tti lam_t];   %时间计加到一块
     lamda_y=[lamdai lamda_y];   %协态累计到一块
              
      if i<N+1           
         inv_lamdai0=[lama(N+1-i+1,1)*(1-p(N+1-i)) lama(N+1-i+1,2)*(1-p(N+1-i))^2 lama(N+1-i+1,3) lama(N+1-i+1,4) lama(N+1-i+1,5)*(1-p(N+1-i))];% 协态跳跃              
       end
 end
tn4(2:end)=cumsum(tn3); %同tn2
 
 %下面把状态和协态进行标准化
 for i=1:N+1;
     
       dt=0.02;  
       stan_t=i-1:dt:i;  %把协态和状态值统一在相同的离散点上

       %下面先把状态标准化
       stan_y1=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),1),stan_t);
       stan_y2=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),2),stan_t);
       stan_y3=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),3),stan_t);
       stan_y4=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),4),stan_t);
       stan_y5=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),5),stan_t);

       %把协态标准化
       stan_lam1=interp1(lam_t(tn4(i)+1:tn4(i+1)),lamda_y(1,tn4(i)+1:tn4(i+1)),stan_t);
       stan_lam2=interp1(lam_t(tn4(i)+1:tn4(i+1)),lamda_y(2,tn4(i)+1:tn4(i+1)),stan_t);
       stan_lam3=interp1(lam_t(tn4(i)+1:tn4(i+1)),lamda_y(3,tn4(i)+1:tn4(i+1)),stan_t);
       stan_lam4=interp1(lam_t(tn4(i)+1:tn4(i+1)),lamda_y(4,tn4(i)+1:tn4(i+1)),stan_t);
       stan_lam5=interp1(lam_t(tn4(i)+1:tn4(i+1)),lamda_y(5,tn4(i)+1:tn4(i+1)),stan_t);
         
        if i<N+1
            grad_p(i)=m*rho1-yb(i,1).*lama(i+1,1)-yb(i,5)*lama(i+1,5)-2*(1-p(i))*yb(i,2)*lama(i+1,2);
         end 
 end
  grad_G=[grad_p]; 
end

%% 非线性约束条件
function [c ceq]=Z_Nostrain(x)
delta=0.1;xi=150;N=5;
u=ODE(x);  
g=Zproj_Var(u);    
c=zeros(5,1);
for i=1:N-1;
    c(i)=g(i+1,1)+norminv(1-delta)*sqrt(g(i+1,2)-g(i+1,1)^2)-xi;
end
ceq=[];
end

%% 投影
function g=Zproj_Var(u)
x0=ones(6,2);
A=[];b=[];Aeq=[];beq=[];VLB=[0 0];VUB=[];
g=fmincon(@(x) Proj_err(x,u),x0,A,b,Aeq,beq,VLB,VUB,@VarConst);
end

function v=Proj_err(x,u)
v=sum(sum((x-u).^2,2));
end

function [c,ceq]=VarConst(x)
c=x(:,1).^2-x(:,2)+.001;   
ceq=[];
end

%% 取出每段终端 N的一阶矩和二阶矩
function u=ODE(h)
T=36;N=5;tau=T/(N+1);
y0=[4 19 3 11 16 ];
for i=1:N+1
      tspan=[i-1:0.1:i];
     [tt yy]=ode45(@Z_state1,tspan,y0,[],tau);
      u(i,:)=[yy(end,1) yy(end,2)];
      if i<N+1
            y0=[yy(end,1)*(1-h(i)) yy(end,2)*(1-h(i))^2  yy(end,3) yy(end,4) yy(end,5)*(1-h(i))];           
      end
end
end
