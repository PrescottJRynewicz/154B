%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prescott Rynewicz, Jordan Robertson, Lukas Kramer
% MAE 154B
% Wing Analysis and Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

% V_n_PDR;

%% Define Wing Shape
M           = 0.02;             %Constants defined for NACA 2412 shape. 
P           = 0.4;
T           = 0.12;
a0          = 0.2969;
a1          = -0.126;
a2          = -0.3516;
a3          = 0.2843;
a4          = -0.1015;

% array for shape
h           =.001;
x_chord     = 0:h:1;
z_camber    = 0:h:1; 
z_thickness = 0:h:1;
for index = 1:length(x_chord)
    if x_chord(index) < P
        z_camber(index) = M/P^2*(2*P*x_chord(index) - x_chord(index)^2);
    else
        z_camber(index) = (M/(1-P)^2)*(1 - 2*P +2*P*x_chord(index) - x_chord(index)^2);
    end
    z_thickness = (T/0.2)*(a0.*x_chord.^.5+a1.*x_chord+a2.*x_chord.^2+a3.*x_chord.^3+a4.*x_chord.^4);
end

upper_surface   = z_camber + z_thickness;
lower_surface   = z_camber - z_thickness;

%% Define Wing Structures
% This section will eventually implement methods to generate different wing
% layouts and positions. 

%Spar Strcuture: Includes the location and area for each spar. Ordered
%clockwise start from top left spar. 
num_spars   = 4; 
spar_pos_1  = 0.2; 
spar_pos_2  = 0.7; 
spar1       = struct('position',[spar_pos_1,0], 'area', 1); 
spar2       = struct('position',[spar_pos_2,0], 'area', 1); 
spar3       = struct('position',[spar_pos_2,0], 'area', 1); 
spar4       = struct('position',[spar_pos_1,0], 'area', 1); 

% Section Struct. Includes the number os stringers, where the section
% start, the web thickness, and arrays of stringer and web objects. 




num_sections    = 4; 
num_stringers   = [5,5,3,3]; 
section1    = struct('num_str', 0,'start_pos', spar1.position(1), 'end_pos', spar2.position(1), 'x_length', 0, 'web_thick', 0, 'stringers', zeros(num_stringers(1),3), 'webs', struct('areas',zeros(num_stringers(1)+1,1)));
section2    = struct('num_str', 0,'start_pos', spar1.position(1), 'end_pos', spar2.position(1), 'x_length', 0, 'web_thick', 0, 'stringers', zeros(num_stringers(2),3), 'webs', struct('areas',zeros(num_stringers(2)+1,1))); 
section3    = struct('num_str', 0,'start_pos', 0                , 'end_pos', spar1.position(1), 'x_length', 0, 'web_thick', 0, 'stringers', zeros(num_stringers(3),3), 'webs', struct('areas',zeros(num_stringers(3)+1,1))); 
section4    = struct('num_str', 0,'start_pos', 0                , 'end_pos', spar1.position(1), 'x_length', 0, 'web_thick', 0, 'stringers', zeros(num_stringers(4),3), 'webs', struct('areas',zeros(num_stringers(4)+1,1)));





% Variables to access stringer arrays. Faster than structures. 
x_pos       = 1; 
z_pos       = 2; 
str_area    = 3; 

%stringer = [x_position, y_position, area]
stringer    = [0,0,0]; 

%web      = start_x_position, end_x_position, thickness
web         = [0,0,0]; 

wing        = struct('spars',[spar1,spar2,spar3,spar4], 'sections', [section1, section2, section3, section4]);


%% Wing Setup: Optimize later to only make updates where changes occur

gap = linspace(0,0,num_sections); 
for index = 1:num_sections
    gap(index)      = (wing.sections(index).end_pos - wing.sections(index).start_pos)/(num_stringers(index)+1); 
end


%% Place all x, z and area values for all sections of wing. 
for index1 = 1:num_sections
    for index2 = 1:num_stringers(index1)
        %% double check index1 = 2
        wing.sections(index1).stringers(index2,x_pos) = wing.sections(index1).start_pos + gap(index1)*index2;
        wing.sections(index1).stringers(index2,z_pos) = get_z(wing.sections(index1).stringers(index2,x_pos),mod(index1,2));
        wing.sections(index1).stringers(index2,str_area) = 1; 
    end
end

%% Place all spar z positions. This will not change, as there will always be 4 spars. 
for index1 = 1:4
    if index1 == 1 || index1 == 2
        wing.spars(index1).position(2) = get_z(wing.spars(index1).position(1),1);
    else
        wing.spars(index1).position(2) = get_z(wing.spars(index1).position(1),0);
    end
end


%% X Centroid

centroid_x  = 0;
centroid_x_area_sum=0;
for index1 = 1:num_sections
    for index2 = 1:num_stringers(index1)
        centroid_x=centroid_x+(wing.sections(index1).stringers(index2,x_pos)...
            *wing.sections(index1).stringers(index2,str_area));
        centroid_x_area_sum=centroid_x_area_sum+wing.sections(index1).stringers(index2,str_area);
    end
end
for index1 = 1:num_spars
    centroid_x = centroid_x + wing.spars(index1).position(x_pos)*wing.spars(index1).area;
    centroid_x_area_sum = centroid_x_area_sum + wing.spars(index1).area;
end


centroid_x = centroid_x/centroid_x_area_sum;

%% Z Centroid
centroid_z=0;
centroid_z_area_sum=0;
for index1 = 1:num_sections
    for index2 = 1:num_stringers(index1)
        centroid_z=centroid_z+(wing.sections(index1).stringers(index2,z_pos)...
            *wing.sections(index1).stringers(index2,str_area));
        centroid_z_area_sum=centroid_z_area_sum+wing.sections(index1).stringers(index2,str_area);
    end
end
for index1 = 1:num_spars
    centroid_z = centroid_z + wing.spars(index1).position(z_pos)*wing.spars(index1).area;
    centroid_z_area_sum = centroid_z_area_sum + wing.spars(index1).area;
end

centroid_z=centroid_z/centroid_z_area_sum;


%% Moment of Inertia
I_x=0;
for index1 = 1:num_sections
    for index2 = 1:num_stringers(index1)
        I_x=I_x+(wing.sections(index1).stringers(index2,str_area)*...
            (wing.sections(index1).stringers(index2,z_pos)-centroid_z)^2);
    end
end

for index1 = 1:num_spars
    I_x = I_x + wing.spars(index1).area*(wing.spars(index1).position(z_pos)-centroid_z)^2; 
end 

I_z=0;
for index1 = 1:num_sections
    for index2 = 1:num_stringers(index1)
        I_z=I_z+(wing.sections(index1).stringers(index2,str_area)*...
            (wing.sections(index1).stringers(index2,x_pos)-centroid_x)^2);
    end
end
for index1 = 1:num_spars
    I_z = I_z + wing.spars(index1).area*(wing.spars(index1).position(x_pos)-centroid_x)^2; 
end 


I_xz=0;
for index1 = 1:num_sections
    for index2 = 1:num_stringers(index1)
        
        I_xz = I_xz + wing.sections(index1).stringers(index2,str_area) *...
            (wing.sections(index1).stringers(index2,z_pos)-centroid_z) *...
            (wing.sections(index1).stringers(index2,x_pos)-centroid_x);
    end
end
for index1 = 1:num_spars
    I_xz = I_xz + wing.spars(index1).area*(wing.spars(index1).position(x_pos)-centroid_x)*...
                  (wing.spars(index1).position(z_pos)-centroid_z); 
end

%% Area cut outs for each stringer
% Double check this for all sections. 
for index1 = 1:num_sections
    for index2 = 1:(num_stringers(index1)+1)
  
        if index2 == 1
            x_plus  = wing.sections(index1).stringers(index2);
            x_minus = wing.sections(index1).start_pos;
            x_i     = wing.sections(index1).start_pos; 
            integral = get_int(x_i,x_plus,mod(index1,2));
        elseif index2 == (num_stringers(index1)+1)
            x_plus  = wing.sections(index1).end_pos;
            x_minus = wing.sections(index1).start_pos;
            x_i     = wing.sections(index1).stringers(index2-1); 
            integral = get_int(x_i,x_plus,mod(index1,2));
        else
            x_plus  = wing.sections(index1).stringers(index2);
            x_minus = wing.sections(index1).start_pos;
            x_i     = wing.sections(index1).stringers(index2-1); 
            integral = get_int(x_i,x_plus,mod(index1,2));
        end
        wing.sections(index1).webs.areas(index2) = abs(get_z(x_i,mod(index1,2))*(x_i-x_minus)/2) - abs((get_z(x_plus,mod(index1,2))*(x_plus-x_minus))/2) + integral(end);     
    end
end

%% Plots

% figure; hold on; axis equal; grid on;
% plot(x_chord,upper_surface,'-')
% plot(x_chord,lower_surface,'-')
% for index = 1:num_sections
%     scatter(wing.sections(index).stringers(:,1),wing.sections(index).stringers(:,2))
% end
% for index = 1:num_spars
%     scatter(wing.spars(index).position(x_pos),wing.spars(index).position(z_pos));
% end
% scatter(centroid_x,centroid_z,'x')
% hold on; xlim([0 1]); ylim([-0.4 0.4]); 
