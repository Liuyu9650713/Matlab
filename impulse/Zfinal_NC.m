%% 编程思路
%  此程序主要采用梯度法求解最优问题，程序架构如下：
%  第一部分  主函数 Zsolve_NC
%                给定优化参数的初始值，
%                给定约束条件，
%                输出优化的参数的结果以及目标函数。
%  第二部分  状态方程 Z_state
%                y(1)->EN
%                y(2)->EN^2
%                y(3)->EC
%                y(4)->EC^2
%                y(5)->ENC
%  第三部分  协态方程  Z_costate
%  第四部分  求目标函数和梯度  Z_objgrade
%              给定参数的取值以及状态的初值
%              通过调用状态方程求出每个时间段的每个时刻的状态值，并在每段末端时刻对 k10 k20 k11进行脉冲作为下一段的初始值
%              根据求得的状态值计算目标函数
%              在求到的状态值的条件下求所对应的协态
%              根据协态求每个想要优化参数的梯度


%%  主函数
function Zfinal_NC
clc ;clear;
N=5;      %农药喷洒次数
x0=zeros(1,N+1);
%x0(1:N)=[0.0027    0.0022    0.0046    0.0060    0.0075];
x0(1:N)=[0.5    0.5    0.5    0.5    0.5];
x0(end)=0.1;

lb=[0*ones(1,N) 0];
ub=[1*ones(1,N) 0.1];
 %A=[1 1 1 1 1 0];
 %b=1;
options = optimset('display','iter');        %迭代的每个值显示出来
options = optimset(options,'GradObj','on');  %用梯度迭代

[x,fval]=fmincon(@Z_objgrade,x0,[],[],[],[],lb,ub,[],options)

end

%% 状态方程
function dy =Z_state(t,y)
alpha=0.0307;lambda=0.0115;eta=0.247e-4;r=0.0026;

dy=zeros(5,1);
dy=[alpha+lambda*y(1)-eta*y(1)*y(3)-eta*y(5);
        alpha+2*lambda*y(2)+lambda*y(1)+eta*y(1)*y(3)-2*eta*y(1)*y(5)-2*eta*y(2)*y(3)+eta*y(5);
        alpha+lambda*y(1)-r*y(3);
        alpha+2*lambda*y(5)+lambda*y(1)-2*r*y(4)+r*y(3);
        alpha+lambda*y(1)+lambda*y(5)+lambda*y(2)-eta*y(1)*y(4)-eta*y(3)*y(5)-r*y(5)];
end

%% 协态方程
function dy=Z_costate(t,y,tt,state1,state2,state3,state4,state5, epsion)
lambda=0.0115;eta=0.247e-4;r=0.0026;delta=0.1;xi=150; d=2;

dy=zeros(5,1);
y1=interp1(tt,state1,t);    %使考虑的点更加全面一些
y2=interp1(tt,state2,t);
y3=interp1(tt,state3,t);
y4=interp1(tt,state4,t);
y5=interp1(tt,state5,t);

dy=(-1)*[2*epsion^(-d)*max(y1+norminv(1-delta)*sqrt(y2)-xi,0)+lambda*y(1)-eta*y(1)*y3+lambda*y(2)+eta*y(2)*y3-2*eta*y(2)*y5+lambda*y(3)+lambda*y(4)+lambda*y(5)-eta*y(5)*y4;
    2*epsion^(-d)*max(y1+norminv(1-delta)*sqrt(y2)-xi,0)*norminv(1-delta)*1/(2*sqrt(y2))+2*lambda*y(2)-2*eta*y(2)*y3+lambda*y(5);
    -eta*y(1)*y1+eta*y(2)*y1-2*eta*y(2)*y2-r*y(3)+r*y(4)-eta*y(5)*y5;
    -2*r*y(4)-eta*y(5)*y1;
    -eta*y(1)-2*eta*y(2)*y1+eta*y(2)+2*lambda*y(4)+lambda*y(5)-eta*y(5)*y3-r*y(5)];
end

%% 目标函数和梯度
function [obj_f,grad_G]=Z_objgrade(x)
rho1=1;rho2=1;A=100;delta=0.1;d=2;beta=1.55;sigma=100;xi=150;
T=168*6;N=5;tau=T/(N+1);

p=zeros(1,N);
p(1:N)=x(1:N);
epsion=x(end);

y0=[10 0.5 10 11 16 ];%各个时间段状态的初始值

% 求目标
ya=zeros(N+1,5);yb=zeros(N+1,5); %把每一段初始值和结束时刻的状态值记录下来 ya:初始  yb:结束
y=[];t=[];                       %空矩阵开始默认是零矩阵，对应求解出来的时间点，状态
tn1=zeros(1,N+1);                %每一段上离散点的个数取长之后构成的行向量，如tn1=[3,5,3,...]第一段上3个离散点，第2段上5个离散点，第3段上3个离散点...
tn2=zeros(1,N+2);                %2*（n+1）+2个时间点上对应的离散点数构成的行向量,如tn2=[0,3,8,11,...]从第2个点开始等于tn1的累加求和
obj_conS=[]; 

 for i=1:N+1;                                         %第一部分：调用状态函数，求解每一点的状态
         tspan=tau*[i-1:0.1:i];
         [tt,yy] = ode45(@Z_state,tspan,y0,[]);   %tt是第i段上离散时间点的列向量，yy状态矩阵，第一列是k10的数量，第二列是k20的数量，形式yy(?,:)，其中？表示第几行，对应离散时刻，：表示第几个状态
         tn1(i)=length(tt);                           %第i段离散时间点取长即点的数目
         
         ya(i,:)=yy(1,:); yb(i,:)=yy(end,:);          %（行向量）把状态在0和1时刻的值记录下来,':'代表第几个状态,'i'代表第i个时间段，1代表0时刻，end代表1时刻，最后累计变成矩阵，第一列k10，第二列k20
         t=[t;tt];                                    %每一次把离散时间点的列向量累计起来，构成列向量t
         y=[y; yy];                                   %把每次的状态累计起来,第一列k10，第二列k20
         
         yyb=max(yy(:,1)+norminv(1-delta)*sqrt(yy(:,2))-xi,0).^2;    %列向量，约束违反函数中的约束 
         obj_conSS=trapz(tt,yyb);                             %约束违反函数积分
         obj_conS=[obj_conS obj_conSS];                           %行向量
         
         if i<N+1
             y0=[yb(i,1)*(1-p(i)) yb(i,2)*(1-p(i))^2 yb(i,3) yb(i,4) yb(i,5)*(1-p(i))];
         end
 end

 tn2(2:end)=cumsum(tn1);                                    %第2个时间点t1到最后一个时间点tN等于每一段的离散点的个数累加求和,如tn2=[0,3,8,11,...，20,23,25]从第2个点开始等于tn1的累加求和

con_viola=sum(obj_conS);                                    %约束违反函数

obj_f=sum(A*rho1*p(1:N))+rho2*yb(end,3)+epsion^(-d)*con_viola+sigma*epsion^beta;

%求梯度
grad_G=zeros(1,N+1);                                        %定义梯度规则
grad_p=zeros(1,N);
grad_epsion=zeros(1,1);
  
 %下面开始调用协态函数求解协态 
lama=zeros(N+1,5);
lamb=zeros(N+1,5);
tn3=zeros(1,N+1);                                           %每一段上离散点的个数取长之后构成的行向量，同tn1

inv_lamdai0=[0 0 rho2 0 0];                                 %协态的在最后一段时间点1上的初始值

 for i=1:N+1;
     
     %i对应第N-i+1个时间段先求第5段，然后第4段，...
     
     ti=t(tn2(N+1-i+1)+1:tn2(N+1-i+2));                     %每一段上的时间点构成的列向量,（对应t里的元素），如果最后一段是两个时间点0.2,0.9，则t1=[0.2 0.9]',...,t5=[0.3 0.4 0.8]',倒求
     yi=y(tn2(N+1-i+1)+1:tn2(N+1-i+2),:);                   %把每一段上每个点的状态记录下来 （对应y里的元素），‘：’代表状态1,状态2,第一列代表S，第2列代表I，倒求
     inv_tspan=[N+1-i+1:-0.1:N+1-i]*tau;

     [inv_tti,inv_lamdai]=ode45(@(t,y)Z_costate(t,y,ti,yi(:,1),yi(:,2),yi(:,3),yi(:,4),yi(:,5),epsion),inv_tspan, inv_lamdai0);
     %inv_tti是第N-i+1段上离散时间点列向量，inv_lamdai是协态的矩阵，第一列是k10的协态，第二列是k20的协态，形式inv_lamdai(?,:)，其中？表示第几个时刻，：表示第几个协态

     tti=fliplr(inv_tti'); lamdai=fliplr(inv_lamdai');              %时间，协态翻转
     lama(N+1-i+1,:)=lamdai(:,1)'; lamb(N+1-i+1,:)=lamdai(:,end)';
     tn3(N+1-i+1)=length(tti);                                      %每一段上对应协态离散化的时间点取长,tn3倒取

      if i<N+1           
         inv_lamdai0=[lama(N+1-i+1,1)*(1-p(N+1-i)) lama(N+1-i+1,2)*(1-p(N+1-i))^2 lama(N+1-i+1,3) lama(N+1-i+1,4) lama(N+1-i+1,5)*(1-p(N+1-i))];% 协态跳跃              
       end
 end
 
 %计算梯度
 for i=1:N+1;        
        if i<N+1
            grad_p(i)=A*rho1-yb(i,1).*lama(i+1,1)-yb(i,5)*lama(i+1,5)-2*(1-p(i))*yb(i,2)*lama(i+1,2);
         end 
 end
  grad_epsion=-d*epsion^(-d-1)*con_viola+sigma*beta*epsion^(beta-1);
  grad_G=[grad_p grad_epsion]; 
end