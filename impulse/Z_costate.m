function dy=Z_costate(t,y,tt,state1,state2,state3,state4,state5, epsion)
a=0.0307;lambda=0.0115;eta=0.247e-4;r=0.0026;rho1=1;rho2=1;m=100;de=0.1;xi=150; d=2;

T=36;N=5;
tau=T/(N+1);

dy=zeros(5,1);
y1=interp1(tt,state1,t);
y2=interp1(tt,state2,t);
y3=interp1(tt,state3,t);
y4=interp1(tt,state4,t);
y5=interp1(tt,state5,t);

dy=(-1)*tau*[2*epsion^(-d)*max(y1+norminv(1-de)*sqrt(y2)-xi,0)+lambda*y(1)-eta*y(1)*y3+lambda*y(2)+eta*y(2)*y3-2*eta*y(2)*y5+lambda*y(3)+lambda*y(4)+lambda*y(5)-eta*y(5)*y4;
    2*epsion^(-d)*max(y1+norminv(1-de)*sqrt(y2-y1^2)-xi,0)*norminv(1-de)*1/(2*sqrt(y2))+2*lambda*y(2)-2*eta*y(2)*y3+lambda*y(5);
    -eta*y(1)*y1+eta*y(2)*y1-2*eta*y(2)*y2-r*y(3)+r*y(4)-eta*y(5)*y5;
    -2*r*y(4)-eta*y(5)*y1;
    -eta*y(1)-2*eta*y(2)*y1+eta*y(2)+2*lambda*y(4)+lambda*y(5)-eta*y(5)*y3-r*y(5)];
end