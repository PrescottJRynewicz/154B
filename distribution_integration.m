Shear Force and Bending Moment for PDR
n=;
w=;
b=;
rho=;
V=;
S=;
Cx=;
y=0:0.001:b;
Fz=(n.*w./b).*((2./pi).*sqrt(1-((2.*y./b).^2))+1);
Fx1=ones(1,round(.8*length(y)));
Fx2=1.25*ones(1,round(.2*length(y)));
Fx=(rho*(V^2)*Cx*S/b)*rot90(horzcat(Fx1,Fx2));
%Bending Moment and Shear Force in X-Direction 
shearx=-cumsum(Fx);
momentx=-cumsum(shearx);
%Bending Moment and Shear Force in z-Direction 
shearz=-cumsum(Fz);
momentz=-cumsum(shearz);