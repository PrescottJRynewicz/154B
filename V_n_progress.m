%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prescott Rynewicz
% 504288967
% MAE 154B Velocity vs. Lift Distrinution for wing loading conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

%% Constants
W                       = 3200;                     % [lb]
b                       = 17*2;                     % [ft]
c                       = 5;                        % [ft]
S                       = 5*17*2;                   % [ft^2]
g                       = 32.2;                     % [ft/s^2]
AR                      = b^2/S;                    % []
rho                     = 0.00238;                  % [imperial] 
mu                      = 1.60E-4;                  % [ft^2/s]
Vc                      = 230;                      % [ft/s]
Vd                      = 270;                      % [ft/s]
Re                      = Vc*c/mu;                  % []
e                       = 0.79;                     % []

%% XFoil Results
xfoil_results = dlmread('naca2412_2.pol');

Cl      = rot90(zeros(length(xfoil_results),1));    % Pre-allocate memory for variables for efficiency.  
alpha   = rot90(zeros(length(xfoil_results),1));     
Cd      = rot90(zeros(length(xfoil_results),1));
Cdp     = rot90(zeros(length(xfoil_results),1));
Cm      = rot90(zeros(length(xfoil_results),1));
min_diff_index  = 1; 



index2 = 1;
for index = 1:length(xfoil_results)
    alpha(index)        = xfoil_results(index,1);   % Pull data from xFoil Results txt file. 
    Cl(index)           = xfoil_results(index,2);
    if (abs(Cl(index))) < abs(Cl(min_diff_index));  % Quick check to find when Cl = 0 for Zero lift Drag condition
        min_diff_index  = index;                    % which is parasite drag. 
    end
    Cd(index)           = xfoil_results(index,3);
    Cdp(index)          = xfoil_results(index,4);
    Cm(index)           = xfoil_results(index,5);
end


%% 2D to 3D Calcs
Cl_alpha                = polyfit(alpha(150:250),...% Get slope of linear portion of Cl. 
                            Cl(150:250),1);
Cl_alpha                = Cl_alpha(1);                  
Cl_alpha_rad            = Cl_alpha*180/pi;
CL_alpha_rad            = Cl_alpha_rad./(1 + ...
                          Cl_alpha_rad./(pi*AR*e)); 
CL                      = CL_alpha_rad/...
                          Cl_alpha_rad*Cl;
CD_i                    = CL.^2./(pi*AR*e); 
CD_p                    = Cd(min_diff_index);       % 0 drag Coefficient: This is the drag coefficient from where CL = 0;
CD                      = (CD_i + CD_p); 
mu_2                    = 2*(W/S)/(rho*g*c*...
                          CL_alpha_rad);
Cn                      = CL.*cosd(alpha) + ...
                          CD.*sind(alpha);
Cx                      = CD.*cosd(alpha) - ...
                          CL.*sind(alpha);               
Cn_alpha_rad            = polyfit(alpha(150:250),Cn(150:250),1); 
Cn_alpha_rad            = Cn_alpha_rad(1)*180/pi;
[Cn_max,max_index]      = max(Cn);
[Cn_min,min_index]      = min(Cn);



%% Dr. Toohey's Code
Cz_max                  = Cn_max;                   %rough estimate for now
Cz_neg                  = Cn_min;                   %rough estimate for now
kg                      = 0.88*mu_2/(5.3+mu_2);
Cza                     = Cn_alpha_rad;             % Lift curve slope


v                       = 0:270;                    % [mph]
vfps                    = v*5280/3600;              % [fps]
n_lift                  = (0.5*rho*S*Cz_max/W).*vfps.^2.*1.25;
n_neg                   =  (0.5*rho*S*Cz_neg/W).*vfps.^2;


lim_load_plus           = 4.4;
lim_load_minus          = -1.76;



g_50_plus               = 1+kg*rho*50*Cza*S/(2*W).*vfps;
g_30_plus               = 1+kg*rho*30*Cza*S/(2*W).*vfps;
g_50_minus              = 1-kg*rho*50*Cza*S/(2*W).*vfps;
g_30_minus              = 1-kg*rho*30*Cza*S/(2*W).*vfps;


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

n_vec                   = [lim_load_plus max(g_50_plus) lim_load_plus min(g_30_minus) min(g_50_minus) -1.76];
v_vec                   = [v(find(n_lift>lim_load_plus,1)) Vc v(end) Vd Vc v(find(n_neg<lim_load_minus,1))]*5280/3600; % [fps]
Cn_vec                  = 2.*n_vec.*W./(rho.*v_vec.^2.*S);

current_diff_index      = 1;
Cx_pos                  = linspace(0,0,6);
for index1 = 1:length(v_vec)
    current_Cn          = Cn_vec(index1);
    current_diff        = 100; 
    current_diff_index  = 1;
    for index2 = 1:length(Cn)
        if (abs((current_Cn - Cn(index2))) < current_diff)
            current_diff_index = index2; 
            current_diff = abs((current_Cn - Cn(index2)));
        end
    end
    Cx_pos(index1) = current_diff_index;
end

Cx_vec                  = (Cx(Cx_pos));

% Fz      = n_vec.*W;          % [lb]
% Fx      = 0.5*rho.*v_vec.^2.*S.*Cx_vec;


%% Code to calc lift distribution, Shear and Bending Moments

y                       =0:0.01:b/2;
length_y                = length(y); 
Fz                      = zeros(length_y,length(Cx_vec)); 
Fx                      = zeros(length_y,length(Cx_vec)); 
shearx                  = zeros(length_y,length(Cx_vec)); 
momentx                 = zeros(length_y,length(Cx_vec));
shearz                  = zeros(length_y,length(Cx_vec)); 
momentz                 = zeros(length_y,length(Cx_vec));




for index=1:length(Cx_vec)
    n                   =n_vec(index);
    V                   =v_vec(index);
    Cx                  =Cx_vec(index);
    Fz(:,index)         =(n.*W./(2.*b)).*((4./pi).*sqrt(1-((2.*y./b).^2))+1);
    Fx1                 =ones(1,round(.8*length(y)));
    Fx2                 =1.25*ones(1,round(.2*length(y)));
    Fx(:,index)         =(rho*(V^2)*Cx*S/b)*horzcat(Fx1,Fx2);
    
    
    %Bending Moment and Shear Force in X-Direction 
    shearx(:,index)     =-cumsum(Fx(:,index));
    momentx(:,index)    =-cumsum(shearx(:,index));
    
    %Bending Moment and Shear Force in z-Direction 
    shearz(:,index)     =-cumsum(Fz(:,index));
    momentz(:,index)    =-cumsum(shearz(:,index));
end



%% Plots 


for index = 1:length(Cx_vec)
    figure; plot(y,Fz(:,index)); title('PHAA F_z Distribution', 'FontSize',16);xlabel('Distace (ft)');ylabel('Force (Lbs)');
    figure; plot(y,shearz(:,index)); title('PHAA V_x Distribution', 'FontSize',16);xlabel('Distace (ft)');ylabel('Force (Lbs)');
    figure; plot(y,momentz(:,index)); title('PHAA Moment_z Distribution', 'FontSize',16);xlabel('Distace (ft)');ylabel('Moment (Lbs*ft)');
    figure; plot(y,Fx(:,index)); title('PHAA F_x Distribution', 'FontSize',16);xlabel('Distace (ft)');ylabel('Force (Lbs)');
    figure; plot(y,shearx(:,index)); title('PHAA V_x Distribution', 'FontSize',16);xlabel('Distace (ft)');ylabel('Force (Lbs)');
    figure; plot(y,momentx(:,index)); title('PHAA M_x Distribution', 'FontSize',16);xlabel('Distace (ft)');ylabel('Moment (Lbs*ft)');
    
end



