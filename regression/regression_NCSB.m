function regression_NCSB
Sbar=200;
k1=0.01;k2=0.01;eta=0.25e-4;alpha=0.03;lambda=0.012;r=0.003;            %程序需要写三层循环，第一层循环是给定每个时刻投放比例和释放比例（若每个时刻给定100组投放比例和释放比例）
mup1=rand(1,100);                                                                               %程序的第二层循环式每次给定一组投放比例和释放比例时模拟100次（即每给定一组释放比例和投放比例使用G算法模拟100次）
mus1=rand(1,100);                                                                                %程序的第三层循环是G算法中的循环（while）
mup2=rand(1,100);
mus2=rand(1,100);
mup3=rand(1,100);
mus3=rand(1,100);
mup4=rand(1,100);
mus4=rand(1,100);
mup5=rand(1,100);
mus5=rand(1,100);
weekhour=[168   336   504   672   840   1008];
S1=[1 1 -1 0;1 1 0 -1;0 0 0 0;0 0 0 0];
S2=[-1 1 0 0;0 0 0 0;-1 1 0 -1;1 -1 -1 0];
S=[S1 S2];
Xa1=[];
xc1=[];
for g=1:100
mpms=[mup1(1,g)  mup2(1,g)  mup3(1,g)  mup4(1,g)  mup5(1,g) 0 mus1(1,g)  mus2(1,g) mus3(1,g)  mus4(1,g)  mus5(1,g) 0];
mp=mpms(1:6);ms=mpms(7:12);
X=[]; % state consisting of n(t), c(t), s(t) and b(t);
X(:,1)=[10;10;0;0];
Obs=[0;10;10;0;0];
TimeNow=0; % time counter, growing until weekhour(1)=168
i=2;
Xb1=[];
for f=1:100
for k=1:6
    while TimeNow<weekhour(k)
        L=size(X,2);
        X_now=X(:,L);
        h=[alpha  lambda*X_now(1)  eta*X_now(1)*X_now(2)  r*X_now(2)  k1*X_now(1)*X_now(3) ...
            k2*X_now(4)  eta*X_now(2)*X_now(4)  eta*X_now(2)*X_now(3)];
        h0=sum(h);
        TimeNow=TimeNow+exprnd(1/h0);
        P=h/h0;
        J=Lookup([1 2 3 4 5 6 7 8],P,1);
        X_New=X(:,i-1)+S(:,J);
        X=[X X_New];
        ObsNew=[TimeNow;X_New];
        Obs=[Obs ObsNew];
        i=i+1;
    end
  L1=length(Obs);
  X(1,L1)=floor((1-mp(k))*X(1,L1-1));
  X(2,L1)=X(2,L1-1)+floor(ms(k)*Sbar);
  X(3,L1)=X(3,L1-1)+floor(ms(k)*Sbar);
  Xb1=[Xb1;X(1,L1)];
end
D1=mean(Xb1(1,:));
D2=mean(Xb1(2,:));
D3=mean(Xb1(3,:));
D4=mean(Xb1(4,:));
D5=mean(Xb1(5,:));
D6=mean(Xb1(6,:));
xc1=[xc1;D1];
xc2=[xc2;D2];
xc3=[xc3;D3];
xc4=[xc4;D4];
xc5=[xc5;D5];
xc6=[xc6;D6];
%L1=length(Obs);
% X(1,L1)=floor((1-mp(k))*X(1,L1-1));
%X(2,L1)=X(2,L1-1)+floor(ms(k)*Sbar);
%X(3,L1)=X(3,L1-1)+floor(ms(k)*Sbar);
end
Xa1=[Xa1;xc1;xc2;xc3;xc4;xc5;xc6];
end
plot(Obs(1,:),X(1,:),'r-')
hold on
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