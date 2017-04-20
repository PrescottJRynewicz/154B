%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prescott Rynewicz
% 504288967
% MAE 154B Velocity vs. Lift Distrinution for wing loading conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

close all;

W = 3200;       % [lb]
b = 17*2;       % [ft]
c = 5;          % [ft]
S = 5*17*2;     % [ft^2]
g = 32.2;       % [ft/s^2]
AR = b^2/S;     % []
rho = 0.00238;  % [imperial] 
mu = 1.60E-4;   % ft^2/s
Vc = 230;       % [ft/s]
Vd = 270;       % [ft/s]
Re = Vc*c/mu;   % []
e = 0.79;       % []

xfoil_results = dlmread('naca2412_2.pol');
xfoil_results = dlmread('dat_data.pol');
Cl      = zeros(length(xfoil_results),1);
alpha   = zeros(length(xfoil_results),1);
Cd      = zeros(length(xfoil_results),1);
Cdp     = zeros(length(xfoil_results),1);
Cm      = zeros(length(xfoil_results),1);
min_diff_index  = 1; 
    


index2 = 1;
for index = 1:length(xfoil_results)
    alpha(index)   = xfoil_results(index,1);
    Cl(index)      = xfoil_results(index,2);
    if (abs(Cl(index))) < abs(Cl(min_diff_index));
        min_diff_index = index;
    end
    Cd(index)      = xfoil_results(index,3);
    Cdp(index)     = xfoil_results(index,4);
    Cm(index)      = xfoil_results(index,5);
end


Cl_alpha        = (Cl(floor(length(Cl)/2)+2) - Cl(floor(length(Cl)/2)-2))/ ...
                    (alpha(floor(length(Cl)/2)+2) - alpha(floor((length(Cl)/2)-2)));  %1/deg
Cl_alpha_rad    = Cl_alpha*180/pi;

CL              = Cl./(1 + Cl./(pi*AR*e)); 
CL_alpha        = (CL(floor(length(CL)/2)+2) - CL(floor(length(CL)/2)-2))/ ...
                    (alpha(floor(length(CL)/2)+2) - alpha((floor(length(CL)/2)-2)));

CL_alpha_rad    = CL_alpha*180/pi;
CD_i            = CL.^2./(pi*AR*e); 
CD_p            = Cd(min_diff_index);  % 0 drag Coefficient: This is the drag coefficient from where CL = 0;
CD              = (CD_i + CD_p); 
mu_2            = 2*(W/S)/(rho*g*c*CL_alpha_rad);


Cn              = CL.*cosd(alpha) + CD.*sind(alpha);
Cx              = CD.*cosd(alpha) - CL.*sind(alpha);

Cn_alpha        = (Cn(floor(length(CL)/2+2)) - Cn(floor(length(CL)/2)-2))/ ...
                    (alpha(floor(length(CL)/2)+2) - alpha((floor(length(CL)/2)-2)));
Cn_alpha_rad    = Cn_alpha*180/pi;


[Cn_max,max_index] = max(Cn);
[Cn_min,min_index] = min(Cn);


Cz_max = Cn_max;            %rough estimate for now
Cz_neg = Cn_min;            %rough estimate for now
kg     = 0.88*mu_2/(5.3+mu_2);
% Cz     = 2*W/(rho*Vc^2*S); 
Cza    = Cn_alpha_rad;               % Lift curve slope


v = 0:270; %mph
vfps = v*5280/3600; %fps
n_lift = (0.5*rho*S*Cz_max/W).*vfps.^2.*1.25;
n_neg =  (0.5*rho*S*Cz_neg/W).*vfps.^2;


lim_load_plus = 4.4;
lim_load_minus = -1.76;



g_50_plus = 1+kg*rho*50*Cza*S/(2*W).*vfps;
g_30_plus = 1+kg*rho*30*Cza*S/(2*W).*vfps;
g_50_minus = 1-kg*rho*50*Cza*S/(2*W).*vfps;
g_30_minus = 1-kg*rho*30*Cza*S/(2*W).*vfps;


figure; grid on; hold on;set(gcf,'color',[1 1 1]);

tempInd = find(n_lift<lim_load_plus); plot(v(tempInd),n_lift(tempInd),'linewidth',2); 
tempInd = find(n_neg>lim_load_minus); plot(v(tempInd),n_neg(tempInd),'linewidth',2);

tempInd = find(n_lift>lim_load_plus,1); plot([v(tempInd) v(end)],[lim_load_plus lim_load_plus],'linewidth',2);

tempInd = find(n_neg<lim_load_minus,1); plot([v(tempInd) Vc],[lim_load_minus lim_load_minus],'linewidth',2);

plot([Vd Vd],[-1 lim_load_plus],'linewidth',2)
plot([Vc Vd],[lim_load_minus -1],'linewidth',2)
plot(v(1:230),g_50_plus(1:230),'--g','linewidth',2)
plot(v(1:270),g_30_plus(1:270),'--g','linewidth',2)
plot(v(1:230),g_50_minus(1:230),'--g','linewidth',2)
plot(v(1:270),g_30_minus(1:270),'--g','linewidth',2)
plot([Vc Vd],[g_50_plus(230) g_30_plus(270)],'--g','linewidth',2)
plot([Vc Vd],[g_50_minus(230) g_30_minus(270)],'--g','linewidth',2)



plot([230 230],[lim_load_minus lim_load_plus],'--r','linewidth',1)
xlabel('V (mph)','fontsize',16,'fontweight','bold');ylabel('n','fontsize',16,'fontweight','bold')
set(gca,'FontSize',16,'fontweight','bold');
ylim([-4 6])



%% Critical Points

n_vec   = [lim_load_plus max(g_50_plus) lim_load_plus min(g_30_minus) min(g_50_minus) -1.76];
v_vec   = [v(find(n_lift>lim_load_plus,1)) Vc v(end) Vd Vc v(find(n_neg<lim_load_minus,1))]*5280/3600; % [fps]
Cn_vec  = 2.*n_vec.*W./(rho.*v_vec.^2.*S);

current_diff_index = 1;
Cx_pos = linspace(0,0,6);
for index1 = 1:length(v_vec)
    current_Cn = Cn_vec(index1);
    current_diff = 100; 
    current_diff_index = 1;
    for index2 = 1:length(Cn)
        if (abs((current_Cn - Cn(index2))) < current_diff)
            current_diff_index = index2; 
            current_diff = abs((current_Cn - Cn(index2)));
        end
    end
    Cx_pos(index1) = current_diff_index;
end

Cx_vec = rot90(Cx(Cx_pos));

Fz      = n_vec.*W;          % [lb]
Fx      = 0.5*rho.*v_vec.^2.*S.*Cx_vec;

