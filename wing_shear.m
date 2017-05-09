%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prescott Rynewicz, Jordan Robertson, Lukas Kramer
% MAE 154B
% Wing Analysis and Optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;

% V_n_PDR;

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
section1    = struct('num_str', 0,'start_pos', spar1.position(1), 'end_pos', spar2.position(1), 'x_length', 0, 'web_thick', 0, 'stringers', zeros(num_stringers(1),3), 'webs', []);
section2    = struct('num_str', 0,'start_pos', spar3.position(1), 'end_pos', spar4.position(1), 'x_length', 0, 'web_thick', 0, 'stringers', zeros(num_stringers(2),3), 'webs', []); 
section3    = struct('num_str', 0,'start_pos', spar4.position(1), 'end_pos', 0                , 'x_length', 0, 'web_thick', 0, 'stringers', zeros(num_stringers(3),3), 'webs', []); 
section4    = struct('num_str', 0,'start_pos', 0                , 'end_pos', spar1.position(1), 'x_length', 0, 'web_thick', 0, 'stringers', zeros(num_stringers(4),3), 'webs', []); 

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


centroid_x=centroid_x/centroid_x_area_sum;

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


%% Plots
xChord = 0:.01:1;
upperSurface = zeros(1,length(xChord));
lowerSurface = zeros(1,length(xChord));

for i=1:length(xChord)
    upperSurface(i) = get_z(xChord(i),1);
    lowerSurface(i) = get_z(xChord(i),0);
end

figure; hold on; axis equal; grid on;
%plot(xChord,z_camber,'-')
plot(xChord,upperSurface,'-')
plot(xChord,lowerSurface,'-')
for index = 1:num_sections
    scatter(wing.sections(index).stringers(:,1),wing.sections(index).stringers(:,2))
end
for index = 1:num_spars
    scatter(wing.spars(index).position(x_pos),wing.spars(index).position(z_pos));
end
scatter(centroid_x,centroid_z,'x')
hold on; xlim([0 1]); ylim([-0.5 0.5]); 
