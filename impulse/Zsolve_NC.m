%% ���˼·
%  �˳�����Ҫ�����ݶȷ�����������⣬����ܹ����£�
%  ��һ����  ������ Zsolve_NC
%                �����Ż������ĳ�ʼֵ��
%                ����Լ��������
%                ����Ż��Ĳ����Ľ���Լ�Ŀ�꺯����
%  �ڶ�����  ״̬���� Z_state
%                y(1)->EN
%                y(2)->EN^2
%                y(3)->EC
%                y(4)->EC^2
%                y(5)->ENC
%  ��������  Э̬����  Z_costate
%  ���Ĳ���  ��Ŀ�꺯�����ݶ�  Z_objgrade
%              ����������ȡֵ�Լ�״̬�ĳ�ֵ
%              ͨ������״̬�������ÿ��ʱ��ε�ÿ��ʱ�̵�״ֵ̬������ÿ��ĩ��ʱ�̶� k10 k20 k11����������Ϊ��һ�εĳ�ʼֵ
%              ������õ�״ֵ̬����Ŀ�꺯��
%              ���󵽵�״ֵ̬��������������Ӧ��Э̬
%              ����Э̬��ÿ����Ҫ�Ż��������ݶ�


%%  ������
function Zsolve_NC
clc ;clear;
N=5;
x0=zeros(1,N+1);
x0(1:N)=[0.2 0.5 0.15 0.2 0.1];
x0(end)=0.1;

lb=[0*ones(1,N) 0];
ub=[1*ones(1,N) 0.1];
% A=[1 1 1 1 1 0];
% b=1;
options = optimset('display','iter');        %������ÿ��ֵ��ʾ����
options = optimset(options,'GradObj','on');  %���ݶȵ���

[x,fval]=fmincon(@Z_objgrade,x0,[],[],[],[],lb,ub,[],options)
end

%% ״̬����
function dy =Z_state(t,y,tau)
alpha=0.0307;lambda=0.0115;eta=0.247e-4;r=0.0026;rho1=1;rho2=1;A=100;delta=0.1;xi=150;

dy=zeros(5,1);
dy=tau*[alpha+lambda*y(1)-eta*y(1)*y(3)-eta*y(5);
        alpha+2*lambda*y(2)+lambda*y(1)+eta*y(1)*y(3)-2*eta*y(1)*y(5)-2*eta*y(2)*y(3)+eta*y(5);
        alpha+lambda*y(1)-r*y(3);
        alpha+2*lambda*y(5)+lambda*y(1)-2*r*y(4)+r*y(3);
        alpha+lambda*y(1)+lambda*y(5)+lambda*y(2)-eta*y(1)*y(4)-eta*y(3)*y(5)-r*y(5)];
end

%% Э̬����
function dy=Z_costate(t,y,tt,state1,state2,state3,state4,state5, epsion)
alpha=0.0307;lambda=0.0115;eta=0.247e-4;r=0.0026;rho1=1;rho2=1;A=100;delta=0.1;xi=150; d=2;

T=168*6;N=5;
tau=T/(N+1);

dy=zeros(5,1);
y1=interp1(tt,state1,t);
y2=interp1(tt,state2,t);
y3=interp1(tt,state3,t);
y4=interp1(tt,state4,t);
y5=interp1(tt,state5,t);

dy=(-1)*tau*[2*epsion^(-d)*max(y1+norminv(1-delta)*sqrt(y2)-xi,0)+lambda*y(1)-eta*y(1)*y3+lambda*y(2)+eta*y(2)*y3-2*eta*y(2)*y5+lambda*y(3)+lambda*y(4)+lambda*y(5)-eta*y(5)*y4;
    2*epsion^(-d)*max(y1+norminv(1-delta)*sqrt(y2)-xi,0)*norminv(1-delta)*1/(2*sqrt(y2))+2*lambda*y(2)-2*eta*y(2)*y3+lambda*y(5);
    -eta*y(1)*y1+eta*y(2)*y1-2*eta*y(2)*y2-r*y(3)+r*y(4)-eta*y(5)*y5;
    -2*r*y(4)-eta*y(5)*y1;
    -eta*y(1)-2*eta*y(2)*y1+eta*y(2)+2*lambda*y(4)+lambda*y(5)-eta*y(5)*y3-r*y(5)];
end

%% Ŀ�꺯�����ݶ�
function [obj_f,grad_G]=Z_objgrade(x)
alpha=0.0307;lambda=0.0115;eta=0.247e-4;r=0.0026;rho1=1;rho2=1;A=100;delta=0.1;d=2;beta=1.55;sigma=100;xi=150;
T=168*6;N=5;tau=T/(N+1);

p=zeros(1,N);
p(1:N)=x(1:N);
epsion=x(end);
p(N+1)=0;

y0=[10 0.5 10 11 16 ];%״̬�ĳ�ʼֵ

% ��Ŀ��
ya=zeros(N+1,5);yb=zeros(N+1,5); %��ÿһ�γ�ʼֵ�ͽ���ʱ�̵�״ֵ̬��¼���� ya:��ʼ  yb:����
y=[];t=[];                       %�վ���ʼĬ��������󣬶�Ӧ��������ʱ��㣬״̬
tn1=zeros(1,N+1);                %ÿһ������ɢ��ĸ���ȡ��֮�󹹳ɵ�����������tn1=[3,5,3,...]��һ����3����ɢ�㣬��2����5����ɢ�㣬��3����3����ɢ��...
tn2=zeros(1,N+2);                %2*n+2��ʱ����϶�Ӧ����ɢ�������ɵ�������,��tn2=[0,3,8,11,...]�ӵ�2���㿪ʼ����tn1���ۼ����
obj_conS=[]; 

 for i=1:N+1;                                         %��һ���֣�����״̬���������ÿһ���״̬
         tspan=[i-1:0.1:i];
         [tt,yy] = ode45(@Z_state,tspan,y0,[],tau);   %tt�ǵ�i������ɢʱ������������yy״̬���󣬵�һ����S���������ڶ�����I����������ʽyy(?,:)�����У���ʾ�ڼ��У���Ӧ��ɢʱ�̣�����ʾ�ڼ���״̬
         tn1(i)=length(tt);                           %ÿһ�α任��[0 1]�ϣ���i����ɢʱ���ȡ���������Ŀ
         
         ya(i,:)=yy(1,:); yb(i,:)=yy(end,:);          %������������״̬��0��1ʱ�̵�ֵ��¼����,':'�����ڼ���״̬,'i'������i��ʱ��Σ�1����0ʱ�̣�end����1ʱ�̣�����ۼƱ�ɾ��󣬵�һ��k10���ڶ���k20
         t=[t;tt];                                    %ÿһ�ΰ���ɢʱ�����������ۼ�����������������t
         y=[y; yy];                                   %��ÿ�ε�״̬�ۼ�����,��һ��k10���ڶ���k20
         
         yyb=max(yy(:,1)+norminv(1-delta)*sqrt(yy(:,2))-xi,0).^2;    %��������Լ��Υ�������е�Լ�� S(t)<xi1��
         obj_conSS=tau*trapz(tt,yyb);                             %Լ��Υ��������S(t)<xi1 �Ļ���
         obj_conS=[obj_conS obj_conSS];                           %������
         
         if i<N+1
             y0=[yb(i,1)*(1-p(i)) yb(i,2)*(1-p(i))^2 yb(i,3) yb(i,4) yb(i,5)*(1-p(i))];
         end
 end

 tn2(2:end)=cumsum(tn1);                                    %��2��ʱ���t1�����һ��ʱ���tN����ÿһ�ε���ɢ��ĸ����ۼ����,��tn2=[0,3,8,11,...��20,23,25]�ӵ�2���㿪ʼ����tn1���ۼ����

con_viola=sum(obj_conS);                                    %Լ��Υ������

obj_f=sum(A*rho1*p(1:N))+rho2*yb(end,3)+epsion^(-d)*con_viola+sigma*epsion^beta;

%���ݶ�
grad_G=zeros(1,N+1);                                        %�����ݶȹ���
grad_p=zeros(1,N);
grad_epsion=zeros(1,1);
  
 %���濪ʼ����Э̬�������Э̬ 
lam_t=[]; lamda_y=[];                                       %�վ���ʼĬ��������󣬶�Ӧ��������ʱ��㣬Э̬
lama=zeros(N+1,5);
lamb=zeros(N+1,5);
tn3=zeros(1,N+1);                                           %ÿһ������ɢ��ĸ���ȡ��֮�󹹳ɵ���������ͬtn1
tn4=zeros(1,N+2);                                           %2*n+2��ʱ����϶�Ӧ����ɢ�������ɵ�������,ͬtn2

inv_lamdai0=[0 0 rho2 0 0];                                 %Э̬�������һ��ʱ���1�ϵĳ�ʼֵ

 for i=1:N+1;
     
     %i��Ӧ��N-i+1��ʱ��������5�Σ�Ȼ���4�Σ�...
     
     ti=t(tn2(N+1-i+1)+1:tn2(N+1-i+2));                     %ÿһ���ϵ�ʱ��㹹�ɵ�������,����Ӧt���Ԫ�أ���������һ��������ʱ���0.2,0.9����t1=[0.2 0.9]',...,t5=[0.3 0.4 0.8]',����
     yi=y(tn2(N+1-i+1)+1:tn2(N+1-i+2),:);                   %��ÿһ����ÿ�����״̬��¼���� ����Ӧy���Ԫ�أ�������������״̬1,״̬2,��һ�д���S����2�д���I������
     inv_tspan=[N+1-i+1:-0.1:N+1-i];

     [inv_tti,inv_lamdai]=ode45(@(t,y)Z_costate(t,y,ti,yi(:,1),yi(:,2),yi(:,3),yi(:,4),yi(:,5),epsion),inv_tspan, inv_lamdai0);
     %inv_tti�ǵ�N-i+1������ɢʱ�����������inv_lamdai��Э̬�ľ��󣬵�һ����S��Э̬���ڶ�����I��Э̬����ʽinv_lamdai(?,:)�����У���ʾ�ڼ���ʱ�̣�����ʾ�ڼ���Э̬

     tti=fliplr(inv_tti'); lamdai=fliplr(inv_lamdai');              %ʱ�䣬Э̬��ת
     lama(N+1-i+1,:)=lamdai(:,1)'; lamb(N+1-i+1,:)=lamdai(:,end)';
     tn3(N+1-i+1)=length(tti);                                      %ÿһ���϶�ӦЭ̬��ɢ����ʱ���ȡ��,tn3��ȡ
     
     lam_t=[tti lam_t];                                             %ʱ��Ƽӵ�һ��
     lamda_y=[lamdai lamda_y];                                      %Э̬�ۼƵ�һ��
              
      if i<N+1           
         inv_lamdai0=[lama(N+1-i+1,1)*(1-p(N+1-i)) lama(N+1-i+1,2)*(1-p(N+1-i))^2 lama(N+1-i+1,3) lama(N+1-i+1,4) lama(N+1-i+1,5)*(1-p(N+1-i))];% Э̬��Ծ              
       end
 end
tn4(2:end)=cumsum(tn3);                                             %ͬtn2
 
 %�����״̬��Э̬���б�׼��
 for i=1:N+1;
     
       dt=0.02;  
       stan_t=i-1:dt:i;                                             %��Э̬��״ֵ̬ͳһ����ͬ����ɢ����

       %�����Ȱ�״̬��׼��
       stan_y1=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),1),stan_t);
       stan_y2=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),2),stan_t);
       stan_y3=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),3),stan_t);
       stan_y4=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),4),stan_t);
       stan_y5=interp1(t(tn2(i)+1:tn2(i+1)),y(tn2(i)+1:tn2(i+1),5),stan_t);

       %��Э̬��׼��
       stan_lam1=interp1(lam_t(tn4(i)+1:tn4(i+1)),lamda_y(1,tn4(i)+1:tn4(i+1)),stan_t);
       stan_lam2=interp1(lam_t(tn4(i)+1:tn4(i+1)),lamda_y(2,tn4(i)+1:tn4(i+1)),stan_t);
       stan_lam3=interp1(lam_t(tn4(i)+1:tn4(i+1)),lamda_y(3,tn4(i)+1:tn4(i+1)),stan_t);
       stan_lam4=interp1(lam_t(tn4(i)+1:tn4(i+1)),lamda_y(4,tn4(i)+1:tn4(i+1)),stan_t);
       stan_lam5=interp1(lam_t(tn4(i)+1:tn4(i+1)),lamda_y(5,tn4(i)+1:tn4(i+1)),stan_t);
         
        if i<N+1
            grad_p(i)=A*rho1-yb(i,1).*lama(i+1,1)-yb(i,5)*lama(i+1,5)-2*(1-p(i))*yb(i,2)*lama(i+1,2);
         end 
 end
  grad_epsion=-d*epsion^(-d-1)*con_viola+sigma*beta*epsion^(beta-1);
  grad_G=[grad_p grad_epsion]; 
end