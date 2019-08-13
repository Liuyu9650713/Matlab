

clear;

  net1=readNet('eee.net');  
  mdyn1=closureDynamics(net1,2,'derMatch','mdyn4_fun');
  [t,mu]=ode23s(@(t,x)mdyn4_fun(x),[0:1:800],mdyn1.Mu.x0);
  mu1=mu(:,1);
  var1=mu(:,3)-mu(:,1).^2;
    
    
% %    plot(t,mu1)
   
   
   %  net1=readNet('cotton_alphid.net');
    [Q,b,c,s,x0]=quadPropensities(net1);
    Q=double(subsParameters(net1,Q)); % r e p l a c e p a r a m e t e r s
    b=double(subsParameters(net1,b)); % by t h e i r n u m e r i c a l v a l u e s
    c=double(subsParameters(net1,c));
    s=double(subsParameters(net1,s));
    nMC=30; % number o f Monte C a r l o r u n s
    Ts=(0:1:800)'; % s i m u l a t i o n t i m e s o f i n t e r e s t
    [X,Xmean,Xstd]=sampledSSA(Q,b,c,s,x0,nMC,Ts);
    
     X1=reshape(X(1,:,:),[nMC length(Ts)]);
     
     subplot(2,1,1);
     plot(Ts,Xmean(1,:));
     subplot(2,1,2);
     plot(Ts,Xstd(1,:).^2);
     
     %plot(Ts,X1,Ts,Xmean(1,:))
     
     %   plot(t,mu1,Ts,Xmean(1,:))
     %   plot(t,var1,Ts,Xstd(1,:).^2)