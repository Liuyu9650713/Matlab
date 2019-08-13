% regression constrain on parameter output B
function TrainData
x=[];C=[];
n=3;
for i=1:n
    ra1=rand(1,5);
    ra2=rand(1,5);
    theta1(i,:)=[ra1/sum(ra1)*rand 0];
    theta2(i,:)=[ra2/sum(ra2)*rand 0];
    theta=[theta1(i,:) theta2(i,:)];% 12 components
    c=ContrF(theta);
    x=[x;theta];
    C=[C;c];
end
input=x% zeros at time 6
output=C% the first column are not to be considered cos no control is conducted
end

function c=ContrF(x)
% x is the input with [mup(1)-mup(6),mus(1)-mus(6)] , the parameters to be controlled
% [c1 c2 c3 c4 c5 c6 c7] are six constraints and moment of C(T) mu1c(T)
% model parameters
Sbar=200;
k1=0.01;k2=0.01;
eta=0.25e-4;alpha=0.03;lambda=0.012;r=0.003;
weekhour=[168   336   504   672   840   1008];
S1=[1 1 -1 0;1 1 0 -1;0 0 0 0;0 0 0 0];S2=[-1 1 0 0;0 0 0 0;-1 1 0 -1;1 -1 -1 0];S=[S1 S2];
mp=x(1:6);ms=x(7:12);
nn=5;
for j=1:nn
    X=[];% state consisting of n(t), c(t), s(t) and b(t);
    X(:,1)= [10;10;0;0];% the first col
    Obs=[0;10;10;0;0];% time and state
    TimeNow=0;% time counter, growing until weekhour(1)=168
    i=2;
    for k=1:6% 6 control times
        while TimeNow<weekhour(k)
            L=size(X,2);
            X_now=X(:,L);% current states
            h=[alpha lambda*X_now(1) eta*X_now(1)*X_now(2) r*X_now(2) k1*X_now(1)*X_now(3) ...
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
        L1=length(Obs);%k=1 the first control time 1
        % size(Obs) say 5,63 so the next is let Obs(:,63)=new impulse Obs=[time and X]
        X(1,L1)=floor((1-mp(k))*X(1,L1-1));% percentage of the ?eld treated with pesticide at time t.
        X(2,L1)=X(2,L1-1)+floor(ms(k)*Sbar);% deterministic amount of steriles added at time t,
        X(3,L1)=X(3,L1-1)+floor(ms(k)*Sbar);
        N(k,j)=X(1,L1-1);%N(k,j)=N(t_k-0,j)
        % the six N(6,j)-N(t_6-0,j)
    end
    CT(j)=X(2,L1-1);
end

for jj=1:6
    c(jj)=log(mean(N(jj,:))+std(N(jj,:))*norminv(0.9));
    % c(2)-c(6) are to be used in fact
end
c(7)=log(mean(CT));
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