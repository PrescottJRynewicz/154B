function z = get_z(x,u)

M = 0.02;
P = 0.4;
T = 0.12;
a0 = 0.2969;
a1 = -0.126;
a2 = -0.3516;
a3 = 0.2843;
a4 = -0.1015;

    if x <P
        z_camber = M/P^2*(2*P*x - x^2);
    else
        z_camber = (M/(1-P)^2)*(1 - 2*P +2*P*x - x^2);
    end

%z_camber = M/P^2*(2*P*x - x^2);
z_thickness = (T/0.2)*(a0*x^.5+a1*x+a2*x^2+a3*x^3+a4*x^4);

if u==1
    z = z_camber + z_thickness;
else
    z = z_camber - z_thickness;
end
