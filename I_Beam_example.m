% Bending/Shear stress example
close all;

length = 20; % in
force = 10000; %lbs

%eye-beam dimensions

max_width = 4; % in
min_width = 1; % in
y_max = 4; % in
center_y = 2; % in


%max bending moment at the root...

M = force*length;

I = min_width*(2*center_y)^3/12 + 2*( max_width*(y_max-center_y)^3/12 + ...
    max_width*(y_max-center_y)*((y_max+center_y)/2)^2);

sigma_max = M * y_max / I;


% solve for shear stress distribution
% V / (I * t) * int(y*da)

% Point 1: evaluated at location just before thickness changes from 4 to 1 in
tempCoeff = force / (I * max_width);
int_y_da = ((y_max+center_y)/2) * max_width*(y_max-center_y);
shear_1 = tempCoeff*int_y_da;

% Point 2: evaluated at location just after thickness changes from 4 to 1 in
tempCoeff = force / (I * min_width);
shear_2 = tempCoeff*int_y_da;


% Point 3: evaluated at center of beam
tempCoeff = force / (I * min_width);
int_y_da = (center_y/2) * min_width*center_y;
shear_3 = shear_2+tempCoeff*int_y_da;

%evaluating continous integral for width of 4..
int_y_da_4 = force / (I * max_width)*4*(y_max^2/2 - (center_y:.1:y_max).^2/2);

%evaluating continous integral for width of 1..
int_y_da_1 = shear_2 + force / (I * min_width)*1*(center_y^2/2 - (0:.1:center_y).^2/2);

figure; grid on; hold on;set(gcf,'color',[1 1 1]);
plot(int_y_da_4,center_y:.1:y_max,'linewidth',2)
plot(int_y_da_1,0:.1:center_y,'linewidth',2)
plot(int_y_da_1,0:-.1:-center_y,'linewidth',2)
plot(int_y_da_4,-center_y:-.1:-y_max,'linewidth',2)
plot([shear_1 shear_2],[center_y center_y],'linewidth',2)
plot([shear_1 shear_2],[-center_y -center_y],'linewidth',2)

plot(shear_1,center_y,'o')
plot(shear_2,center_y,'o')
plot(shear_3,0,'o')
plot(shear_2,-center_y,'o')
plot(shear_1,-center_y,'o')

xlabel('shear stress (lb/in^2)','fontsize',16,'fontweight','bold');ylabel('Distance from Center (in)','fontsize',16,'fontweight','bold')
set(gca,'FontSize',16,'fontweight','bold');


figure; grid on; hold on;set(gcf,'color',[1 1 1]);
plot([0 4 4 2.5 2.5 4 4 0 0 1.5 1.5 0 0],[4 4 2 2 -2 -2 -4 -4 -2 -2 2 2 4],'linewidth',2)

