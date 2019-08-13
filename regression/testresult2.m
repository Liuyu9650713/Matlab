function testresult2
% ODE test of control
N=6;
y0=[10 0 10 0 0];
%h=[0.99 0.99 0.99 0.99 0.99 0];
ra1=rand(1,5);
h=[ra1/sum(ra1)*rand 0];
y1=[];
y2=[];
t1=[];
for i=1:N
    tspan=[168*(i-1) 168*i];
    [tt yy]=ode45(@state_NC,tspan,y0,[]);
    t1=[t1;tt];
    y1=[y1;yy(:,1)];
    %y2=[y2;yy(:,2)];
    y0=[yy(end,1)*(1-h(i)) yy(end,2) yy(end,3)*(1-h(i))^2 yy(end,4)*(1-h(i)) yy(end,5)];
end
hold on
plot(t1,y1,'k-')
box on
% 1.EN;2.EC;3.EN^2;4.ENC;5.EC^2
% G simulation of control 
eta=0.25e-4;alpha=0.03;lambda=0.012;r=0.003;
weekhour=[168   336   504   672   840   1008];%observinbg times
S=[1 1 -1 0;1 1 0 -1]; 
mp=h;
X=[];% state consisting of n(t), c(t), s(t) and b(t);
X(:,1)= [10;10];% the first col
Obs=[0;10;10];% time and state
TimeNow=0;% time counter, growing until weekhour(1)=168
i=2;
for k=1:6
    while TimeNow<weekhour(k)
        L=size(X,2);
        X_now=X(:,L);% current states
        h=[alpha    lambda*X_now(1)   eta*X_now(1)*X_now(2)    r*X_now(2)];% computing h(x,ci) at x
        h0=sum(h);
        TimeNow=TimeNow+exprnd(1/h0);% time to response
        P=h/h0; J=Lookup([1 2 3 4],P,1);% next response
        X_New=X(:,i-1)+S(:,J);% new states
        X=[X X_New];%adding by states
        ObsNew=[TimeNow;X_New];% new time and two states
        Obs=[Obs ObsNew];% adding by new time and two states
        i=i+1;
    end
    L1=length(Obs);
    X(1,L1)=floor((1-mp(k))*X(1,L1-1));% percentage of the field treated with pesticide at time t.
end
hold on
%plot(Obs(1,:),X(1,:),'r-')
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


function dx =state_NC(t,x) 
%alpha=0.03;lambda=0.012;eta=0.000025;r=0.003; 
% 1.EN;2.EC;3.EN^2;4.ENC;5.EC^2
% x(1)表示N的一阶矩，x(2)表示C的一阶矩，x(3)表示N的二阶矩，
% X(4)表示NC的混合矩，X(5)表示C的二阶矩，
alpha=0.0307;lambda=0.0115;eta=0.0000247;r=0.0026; 
dx=zeros(5,1);
  dx(1,1) = alpha - eta*x(4) + lambda*x(1);
  dx(2,1) = alpha + lambda*x(1) - r*x(2);
  dx(3,1) = alpha + eta*x(4) + 2*lambda*x(3) + x(1)*(2*alpha + lambda) - 2*eta*(2*x(1)*x(4) + x(2)*x(3) - 2*x(1)^2*x(2));
  dx(4,1) = alpha + x(1)*(alpha + lambda) + alpha*x(2) + lambda*x(3) - eta*(x(1)*x(5) + 2*x(2)*x(4) - 2*x(1)*x(2)^2) + x(4)*(lambda - r);
  dx(5,1) = alpha + lambda*x(1) + 2*lambda*x(4) - 2*r*x(5) + x(2)*(2*alpha + r);
end


