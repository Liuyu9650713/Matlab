function GandODEcontrast
% ODE test of control
N=6;
y0=[10 0 10 0 0];
h=[0 0 0 0 0 0];% possible control amounts
y1=[];t1=[];
%Conthours=[168   336   504   672   840   1008];
for i=1:N
    tspan=[168*(i-1) 168*i];
    [tt yy]=ode45(@state_NC,tspan,y0,[]);
    t1=[t1;tt];
    y1=[y1;yy(:,1)];
    y0=[yy(end,1)*(1-h(i)) yy(end,2) yy(end,3) yy(end,4) yy(end,5)];
end
hold on
plot(t1,y1(:,1))
box on
% G simulation of control 
Sbar=200;% A=100;

k1=0.01;k2=0.01;
eta=0.247e-4;alpha=0.0307;lambda=0.0115;r=0.0026;


%k1=0.01;k2=0.01;eta=0.25e-4;alpha=0.03;lambda=0.012;r=0.003;
weekhour=[168   336   504   672   840   1008];%observinbg times
S1=[1 1 -1 0;1 1 0 -1;0 0 0 0;0 0 0 0];
S2=[-1 1 0 0;0 0 0 0;-1 1 0 -1;1 -1 -1 0];
S=[S1 S2];
%mpms=[0    0.0000    0.2941         0         0         0    0.1162         0         0         0         0         0];
mpms=[ 0    0.4001    0.1359    0.6982    0.5658 0  ...
    0.3548 0 0 0 0 0];


mp=mpms(1:6);ms=mpms(7:12);
X=[];% state consisting of n(t), c(t), s(t) and b(t);
X(:,1)= [10;10;0;0];% the first col
Obs=[0;10;10;0;0];% time and state
TimeNow=0;% time counter, growing until weekhour(1)=168
i=2;
for k=1:6
    while TimeNow<weekhour(k)
        L=size(X,2);
        X_now=X(:,L);% current states
        h=[alpha    lambda*X_now(1)   eta*X_now(1)*X_now(2)    r*X_now(2)     k1*X_now(1)*X_now(3) ...
            k2*X_now(4)  eta*X_now(2)*X_now(4)  eta*X_now(2)*X_now(3)];% computing h(x,ci) at x
        h0=sum(h);
        TimeNow=TimeNow+exprnd(1/h0);% time to response
        P=h/h0; J=Lookup([1 2 3 4 5 6 7 8],P,1);% next response
        X_New=X(:,i-1)+S(:,J);% new states
        X=[X X_New];%adding by states
        ObsNew=[TimeNow;X_New];% new time and two states
        Obs=[Obs ObsNew];% adding by new time and two states
        i=i+1;
    end
    L1=length(Obs);
    X(1,L1)=floor((1-mp(k))*X(1,L1-1));% percentage of the ?eld treated with pesticide at time t.
    X(2,L1)=X(2,L1-1)+floor(ms(k)*Sbar);% deterministic amount of steriles added at time t,
    X(3,L1)=X(3,L1-1)+floor(ms(k)*Sbar);
end
hold on
plot(Obs(1,:),X(1,:),'r-')
plot(Obs(1,:),X(2,:),'k-')
% xlabel('Time(hrs)','fontsize',10)
% ylabel('E[N(t)]','fontsize',10)
box on
end


function samp=Lookup(Y,P,n)
% P may have zeros in its components
P=P(:);Y=Y(:);PY=[P Y];PY=PY(P~=0,:);
P=PY(:,1);Y=PY(:,2);
% then P with no zeros in its components
F=cumsum(P);
for i=1:n
    k=1;
    u=rand;
    while F(k)<u
        k=k+1;
    end
    x(i)=k;
end
samp=Y(x,:);
end

function dy =state_NC(t,y) %考虑只喷洒农药的矩方程
alpha=0.0307;lambda=0.0115;eta=0.0000247;r=0.0026; 
% y(1)表示N的一阶矩，y(2)表示N的二阶矩，y(3)表示C的一阶矩，y(4)表示C的二阶矩，y(5)表示NC的混合矩
dy=zeros(5,1);
dy=[alpha+lambda*y(1)-eta*y(1)*y(3)-eta*y(5);
    alpha+2*lambda*y(2)+lambda*y(1)+eta*y(1)*y(3)-2*eta*y(1)*y(5)-2*eta*y(2)*y(3);
    alpha+lambda*y(1)-r*y(3);
    alpha+2*lambda*y(5)+lambda*y(1)-2*r*y(4)+r*y(3);
    alpha+lambda*y(1)+lambda*y(5)+lambda*y(2)-eta*y(1)*y(4)-eta*y(3)*y(5)-eta*(y(1)*y(4)+2*y(3)*y(5)+y(1)*y(3)^2)-r*y(5)];
end