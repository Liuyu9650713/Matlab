%the trace of N and C without  control 
function TraceWithoutControl
alpha=0.03;lambda=0.012;eta=0.25*10^(-4);r=0.003;
weekhour=[168   336   504   672   840  1000];%observinbg times
%z=lognrnd(-0.0244,0.2209);
    for j=1:2% 100 plants used for obs
        X=[];% state consisting of n(t) and c(t); two long longs
        X(:,1)=[10 10]';% the first col
        Obs=[0;10;10];% time and state
        TimeNow=0;% time counter, growing until weekhour(5)=840
        S=[1 1 -1 0;1 1 0 -1];%stoichiometry matrix
        i=2;
        while TimeNow<weekhour(end)
            L=size(X,2);
            X_now=X(:,L);% current states
            z=lognrnd(-0.0244,0.2209);
            h=[alpha lambda*X_now(1)*z eta*X_now(1)*X_now(2) r*X_now(2)];% computing h(x,ci) at x
            h0=sum(h);
            TimeNow=TimeNow+exprnd(1/h0);% time to response
            P=h/h0;
            J=Lookup([1 2 3 4],P,1);% next response
            X_New=X(:,i-1)+S(:,J);% new states
            X=[X X_New];%adding by states
            ObsNew=[TimeNow;X_New];% new time and two states
            Obs=[Obs ObsNew];% adding by new time and two states
            i=i+1;
        end
    end
plot(Obs(1,:),X(1,:))
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