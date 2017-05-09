function z = get_int(xi,xf,u)

M = 0.02;
P = 0.4;
T = 0.12;
a0 = 0.2969;
a1 = -0.126;
a2 = -0.3516;
a3 = 0.2843;
a4 = -0.1015;


%evaluate the integral of camber line, depending on xi and xf related to P

if xf <P
    intCamb =  M/P^2*(2*P*xf^2/2 - xf^3/3) - M/P^2*(2*P*xi^2/2 - xi^3/3);   
elseif xi<P
    intCamb = (M/(1-P)^2)*((1 - 2*P)*xf +2*P*xf^2/2 - xf^3/3) - (M/(1-P)^2)*((1 - 2*P)*P +2*P*P^2/2 - P^3/3);
    intCamb = intCamb + M/P^2*(2*P*P^2/2 - P^3/3) - M/P^2*(2*P*xi^2/2 - xi^3/3);
else
    intCamb = (M/(1-P)^2)*((1 - 2*P)*xf +2*P*xf^2/2 - xf^3/3) - (M/(1-P)^2)*((1 - 2*P)*xi +2*P*xi^2/2 - xi^3/3);
end

% do integral on thickness line
%z_thickness = (T/0.2)*(a0*x^.5+a1*x+a2*x^2+a3*x^3+a4*x^4);

intThickness = (T/0.2)*(a0*xf^1.5/1.5 + a1*xf^2/2 + a2*xf^3/3 + a3*xf^4/4 +a4*xf^5/5);
intThickness = intThickness - (T/0.2)*(a0*xi^1.5/1.5 + a1*xi^2/2 + a2*xi^3/3 + a3*xi^4/4 +a4*xi^5/5);

% combine both integral results to get total integral
if u == 1
   z = intCamb + intThickness; 
else
   z = abs(intCamb - intThickness);
end