close all;
force = 8000; % lbs
stringer_A = 0.5; % in^2
thickness = 0.04; % in

top_stringers_y = 6; % in
middle_stringers_y = 2; % in

I = 2*stringer_A*top_stringers_y^2 + 2*stringer_A*middle_stringers_y^2;

% solve for shear stress distribution.  this calc ignores the thickness of
% the web between teh stringers (assumes bending taken by stringers)
% V / (I * t) * int(y*da)

shear_top_web = force / (I*thickness) * top_stringers_y * stringer_A;
shear_middle_web = shear_top_web + (force / (I*thickness) * middle_stringers_y * stringer_A);


figure; grid on; hold on;set(gcf,'color',[1 1 1]);


plot([shear_top_web shear_top_web],[middle_stringers_y top_stringers_y],'linewidth',2);
plot([shear_middle_web shear_middle_web],[-middle_stringers_y middle_stringers_y],'linewidth',2);
plot([shear_top_web shear_top_web],[-middle_stringers_y -top_stringers_y],'linewidth',2);

plot([0 shear_top_web],[top_stringers_y top_stringers_y],'linewidth',2);
plot([0 shear_top_web],[-top_stringers_y -top_stringers_y],'linewidth',2);
plot([shear_middle_web shear_top_web],[middle_stringers_y middle_stringers_y],'linewidth',2);
plot([shear_middle_web shear_top_web],[-middle_stringers_y -middle_stringers_y],'linewidth',2);
xlabel('shear stress (lb/in^2)','fontsize',16,'fontweight','bold');ylabel('Distance from Center (in)','fontsize',16,'fontweight','bold')
set(gca,'FontSize',16,'fontweight','bold');

%Alternate approach.. compute change in bending stress at each stringer to
%find the change in shear load

%at top stringer
d_sigma = force * top_stringers_y / I;  %(lbs/in^2)
d_force_top = d_sigma * stringer_A;

%at middle stringer..
d_sigma = force * middle_stringers_y / I;  %(lbs/in^2)
d_force_middle = d_force_top + d_sigma*stringer_A;

%check if load balances
check_load = 2*d_force_top*4 + d_force_middle*4;


