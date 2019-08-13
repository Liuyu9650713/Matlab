clc;clear;
a=0.0307;lambda=0.0115;eta=0.247e-4;r=0.0026;rho1=1;rho2=1;m=100;de=0.1;xi=150; d=2;b=1.55;s=100;

T=36;N=5;tau=T/(N+1);
p=[0 0.25 0.3 0.2 0.1];%时间段theta(i)
dt=0.001;p1=0:dt:0.15;
epsion=0.1;

ya=zeros(N,5);yb=zeros(N,5);%把每段的初始和结束时刻都记录下来，ya记录初始状态，yb代表结束状态

tn1=zeros(1,N+1);
tn2=zeros(1,N+2);

obj_f1=[];
grad_p1=[];
obj_conS=[];

for j=1:0.15/dt+1;
    y0=[10 19 10 11 16];
    y=[];t=[];
    tn1=zeros(1,N+1);
    tn2=zeros(1,N+2);
    
    for i=1:N+1;%第一部分：调用状态函数，求解每一点的状态
         tspan=[i-1:0.1:i];
         [tt,yy] = ode45(@Z_state,tspan,y0,[],tau);%tt是第i段上离散时间点的列向量，yy状态矩阵，第一列是S的数量，第二列是I的数量，形式yy(?,:)，其中？表示第几行，对应离散时刻，：表示第几个状态
         tn1(i)=length(tt);%每一段变换到[0 1]上，第i段离散时间点取长即点的数目
         
         tn1(i)=length(tt);%每一段变换到[0 1]上，第i段离散时间点取长即点的数目
         u=yy(:,1:2);
         g=Proj_Var(u,length(tt));
         yy(:,1)=g(:,1);
         yy(:,2)=g(:,2);
         
         ya(i,:)=yy(1,:); yb(i,:)=yy(end,:); %（行向量）把状态在0和1时刻的值记录下来,':'代表第几个状态,'i'代表第i个时间段，1代表0时刻，end代表1时刻，最后累计变成矩阵，第一列S，第二列I
         t=[t;tt]; %每一次把离散时间点的列向量累计起来，构成列向量t
         y=[y; yy]; %  y is a length(tt)*2 vector; 把每次的状态累计起来,第一列S，第二列I
          
         yyb=max(yy(:,1)+norminv(1-de)*sqrt(yy(:,2))-xi,0).^2;%列向量，约束违反函数中的约束 S(t)<xi1，
         obj_conSS=tau*trapz(tt,yyb);%数，约束违反函数中S(t)<xi1 的积分
         obj_conS=[obj_conS obj_conSS];%行向量
         
         if i<N+1
             switch i
                 case 1
                    y0=[yb(i,1)*(1-p1(j)) yb(i,2)*(1-p1(j))^2 yb(i,3) yb(i,4) yb(i,5)*(1-p1(j))];
                 case 2
                     y0=[yb(i,1)*(1-p(i)) yb(i,2)*(1-p(i))^2 yb(i,3) yb(i,4) yb(i,5)*(1-p(i))];
                 case 3
                     y0=[yb(i,1)*(1-p(i)) yb(i,2)*(1-p(i))^2 yb(i,3) yb(i,4) yb(i,5)*(1-p(i))];
                 case 4
                     y0=[yb(i,1)*(1-p(i)) yb(i,2)*(1-p(i))^2 yb(i,3) yb(i,4) yb(i,5)*(1-p(i))];
                 case 5
                     y0=[yb(i,1)*(1-p(i)) yb(i,2)*(1-p(i))^2 yb(i,3) yb(i,4) yb(i,5)*(1-p(i))];
             end
         end
    end
 

 tn2(2:end)=cumsum(tn1); %第2个时间点t1到最后一个时间点tN等于每一段的离散点的个数累加求和,如tn2=[0,3,8,11,...，20,23,25]从第2个点开始等于tn1的累加求和

con_viola=sum(obj_conS);
obj_f=sum(m*rho1*p(2:N))+m*rho1*p1(j)+rho2*yb(end,3)+epsion^(-d)*con_viola+s*epsion^b;
obj_f1=[obj_f1 obj_f];

lam_t=[]; lamda_y=[];%空矩阵开始默认是零矩阵，对应求解出来的时间点，协态
lama=zeros(N+1,5);
lamb=zeros(N+1,5);
tn3=zeros(1,N+1); %每一段上离散点的个数取长之后构成的行向量，同tn1
tn4=zeros(1,N+2);

inv_lamdai0=[0 0 rho2 0 0];

 for i=1:N+1;
     
     %i对应第N-i+1个时间段先求第5段，然后第4段，...
     
     ti=t(tn2(N+1-i+1)+1:tn2(N+1-i+2));%每一段上的时间点构成的列向量,（对应t里的元素），如果最后一段是两个时间点0.2,0.9，则t1=[0.2 0.9]',...,t5=[0.3 0.4 0.8]',倒求
     yi=y(tn2(N+1-i+1)+1:tn2(N+1-i+2),:);%把每一段上每个点的状态记录下来 （对应y里的元素），‘：’代表状态1,状态2,第一列代表S，第2列代表I，倒求
     inv_tspan=[N+1-i+1:-0.1:N+1-i];

     [inv_tti,inv_lamdai]=ode45(@(t,y)Z_costate(t,y,ti,yi(:,1),yi(:,2),yi(:,3),yi(:,4),yi(:,5),epsion),inv_tspan, inv_lamdai0);
     %inv_tti是第N-i+1段上离散时间点列向量，inv_lamdai是协态的矩阵，第一列是S的协态，第二列是I的协态，形式inv_lamdai(?,:)，其中？表示第几个时刻，：表示第几个协态

     tti=fliplr(inv_tti'); lamdai=fliplr(inv_lamdai'); %时间，协态翻转
     lama(N+1-i+1,:)=lamdai(:,1)'; lamb(N+1-i+1,:)=lamdai(:,end)';
     tn3(N+1-i+1)=length(tti);%每一段上对应协态离散化的时间点取长,tn3倒取
     
     lam_t=[tti lam_t];   %时间计加到一块
     lamda_y=[lamdai lamda_y];   %协态累计到一块
              
      if i<N+1                                  
          switch N+1-i
                 case 1
                    inv_lamdai0=[lama(N+1-i+1,1)*(1-p1(j)) lama(N+1-i+1,2)*(1-p1(j))^2 lama(N+1-i+1,3) lama(N+1-i+1,4) lama(N+1-i+1,5)*(1-p1(j))];% 协态跳跃
                 case 2
                     inv_lamdai0=[lama(N+1-i+1,1)*(1-p(N+1-i)) lama(N+1-i+1,2)*(1-p(N+1-i))^2 lama(N+1-i+1,3) lama(N+1-i+1,4) lama(N+1-i+1,5)*(1-p(N+1-i))];% 协态跳跃
                 case 3
                     inv_lamdai0=[lama(N+1-i+1,1)*(1-p(N+1-i)) lama(N+1-i+1,2)*(1-p(N+1-i))^2 lama(N+1-i+1,3) lama(N+1-i+1,4) lama(N+1-i+1,5)*(1-p(N+1-i))];% 协态跳跃
                 case 4
                     inv_lamdai0=[lama(N+1-i+1,1)*(1-p(N+1-i)) lama(N+1-i+1,2)*(1-p(N+1-i))^2 lama(N+1-i+1,3) lama(N+1-i+1,4) lama(N+1-i+1,5)*(1-p(N+1-i))];% 协态跳跃
                 case 5
                    inv_lamdai0=[lama(N+1-i+1,1)*(1-p(N+1-i)) lama(N+1-i+1,2)*(1-p(N+1-i))^2 lama(N+1-i+1,3) lama(N+1-i+1,4) lama(N+1-i+1,5)*(1-p(N+1-i))];% 协态跳跃
             end
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
 end         
        grad_p=m*rho1-yb(1,1).*lama(2,1)-yb(1,5)*lama(2,5)-2*(1-p1(j))*yb(1,2)*lama(2,2);
        grad_p1=[grad_p1 grad_p];
end

dt=0.001;
obj_fgrad=diff(obj_f1)/dt;%行向量，差导，25维
                
figure(1)
plot(p1,obj_f1,'b','LineWidth',2)
hold on
figure(2)
plot(p1,grad_p1,'ko',p1(1:0.15/dt),obj_fgrad,'r+','LineWidth',1)%红色的是差分
hold on
