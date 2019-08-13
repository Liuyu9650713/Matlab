
a=0.0307;lambda=0.0115;eta=0.247e-4;r=0.0026;rho1=1;rho2=1;m=100;de=0.1;xi=150; d=2;b=1.55;s=100;
N=5; T=168*6;

%��һ�飺nonoptimal
y0=[10 5 10 11 16 ];    


x1(1:N)=[0.5 0.5 0.5 0.5 0.5];


tau=T/(N+1); p=zeros(1,N);
p(1:N)=x1(1:N);
Tspan1=[0:0.01:1];

ya=zeros(N+1,5);yb=zeros(N+1,5);
y1=[];t1=[];
tn1=zeros(1,N+1);
tn2=zeros(1,N+2);

for i=1:N+1; %��һ�飺�̶�����
    
    [tt1,yy1]=ode45(@Z_state,Tspan1,y0,[],tau); %solve state, tt is the column vector, yy is the matrix(column:S,I,P,V)
    tn1(i)=length(tt1); %the number of the points
    ya(i,:)=yy1(1,:); yb(i,:)=yy1(end,:);  %��״̬�ڳ�ʼ�ͽ���ʱ�̵�ֵ��¼������':'����ڼ���״̬,'i'����ڼ���ʱ�����
    
    tsw1=[0 tau 2*tau 3*tau 4*tau 5*tau];
    t1=[t1; tsw1(i)+tt1*tau];%��ԭ��ԭ����ʱ���ϣ���ʱ���ۼ�������������
    y1=[y1; yy1];%��״̬�ۼ���������һ��S���ڶ���I
    
    
    if i<N+1
        y0=[yb(i,1)*(1-p(i)) yb(i,2)*(1-p(i))^2 yb(i,3) yb(i,4) yb(i,5)*(1-p(i))]; %״̬����Ծ

    end
       
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%�ڶ��飺�Ż�����

y0=[10 5 10 11 16 ];    %״̬N,S,R,SR,M�ĳ�ʼֵ

x2(1:N)=[ 0.6882    0.7159    0.6470    0.6436    0.4617];
%x2(1:N)=[0.6948    0.8018    0.9584    0.6691    0.7376];

tau=T/(N+1); p=zeros(1,N);
p(1:N)=x2(1:N);

epsion=x2(end);

ya2=zeros(N+1,5);yb2=zeros(N+1,5);
y2=[];t2=[];
tn1=zeros(1,N+1);
tn2=zeros(1,N+2);  
Tspan2=[0:0.01:1];%ʱ��任��[i-1,i]
for i=1:N+1;
    
    %i=1��Ӧ��1��
    
    [tt2,yy2]=ode45(@Z_state,Tspan2,y0,[],tau); %solve state, tt is the column vector, yy is the matrix(column:S,I,P,V)
    tn2(i)=length(tt2); %the number of the points
    
    ya2(i,:)=yy2(1,:); yb2(i,:)=yy2(end,:);  %��״̬�ڳ�ʼ�ͽ���ʱ�̵�ֵ��¼������':'����ڼ���״̬,'i'����ڼ���ʱ�����
    
    tsw2=[0 tau 2*tau 3*tau 4*tau 5*tau];
    t2=[t2; tsw2(i)+tt2*tau];%��ԭ��ԭ����ʱ���ϣ���ʱ���ۼ�������������
    y2=[y2; yy2];%��״̬�ۼ���������һ��S���ڶ���I
    
     if i<N+1
           y0=[yb2(i,1)*(1-p(i)) yb2(i,2)*(1-p(i))^2  yb2(i,3) yb2(i,4) yb2(i,5)*(1-p(i))];%��Ӷ�Ԫ������״̬����Ծ           
     end  
end


% figure(1)
% plot(t1,y1(:,1),'b:',t2,y2(:,1),'r','LineWidth',1.5)
% xlabel('Time(hours)','FontSize',12);
% ylabel('k10','FontSize',12);
% legend('Before optimization','After optimization');
% set(gca,'linewidth',2);
figure(1)
plot(t1,y1(:,3),'b:',t2,y2(:,3),'r','LineWidth',1.5)
xlabel('Time(hours)','FontSize',12);
ylabel('k01','FontSize',12);
legend('Before optimization','After optimization');
set(gca,'linewidth',2);
 hold on

% figure(2)
% plot(t1,y1(:,2),'b:',t2,y2(:,2),'r','LineWidth',1)
% xlabel('Time(hours)','FontSize',12);
% ylabel('k20','FontSize',12);
% 
% hold on
% 
% figure(3)
% plot(t1,y1(:,3),'b:',t2,y2(:,3),'r','LineWidth',1.5)
% xlabel('Time(hours)','FontSize',12);
% ylabel('k01','FontSize',12);
% hold on
% 
% figure(4)
% plot(t1,y1(:,4),'b:',t2,y2(:,4),'r','LineWidth',1.5)
% xlabel('Time(hours)','FontSize',12);
% ylabel('k02','FontSize',12);
% 
% hold on
% 
% figure(5)
% plot(t1,y1(:,5),'b:',t2,y2(:,5),'r','LineWidth',1.5)
% xlabel('Time(hours)','FontSize',12);
% ylabel('k11','FontSize',12);
% legend('Before optimization','After optimization');
% set(gca,'linewidth',2);
% hold on
