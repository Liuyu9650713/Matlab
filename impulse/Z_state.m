function dy =Z_state(t,y,tau)
a=0.0307;lambda=0.0115;eta=0.247e-4;r=0.0026;rho1=1;rho2=1;m=100;de=0.1;xi=150;

dy=zeros(5,1);
dy=tau*[a+lambda*y(1)-eta*y(1)*y(3)-eta*y(5);
        a+2*lambda*y(2)+lambda*y(1)+eta*y(1)*y(3)-2*eta*y(1)*y(5)-2*eta*y(2)*y(3)+eta*y(5);
        a+lambda*y(1)-r*y(3);
        a+2*lambda*y(5)+lambda*y(1)-2*r*y(4)+r*y(3);
        a+lambda*y(1)+lambda*y(5)+lambda*y(2)-eta*y(1)*y(4)-eta*y(3)*y(5)-r*y(5)];
end