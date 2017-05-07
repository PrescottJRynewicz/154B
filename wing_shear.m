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
spar1       = struct('position',[0.20,0], 'area', 0); 
spar2       = struct('position',[0.75,0], 'area', 0); 
spar3       = struct('position',[0.75,0], 'area', 0); 
spar4       = struct('position',[0.20,0], 'area', 0); 

% Section Struct. Includes the number os stringers, where the section
% start, the web thickness, and arrays of stringer and web objects. 
num_sections    = 4; 
num_stringers   = [10,2,2,2]; 
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

for index1 = 1:num_sections
    for index2 = 1:num_stringers(index1)
        %% double check index1 = 2
        wing.sections(index1).stringers(index2,x_pos) = wing.sections(index1).start_pos + gap(index1)*index2;
        wing.sections(index1).stringers(index2,z_pos) = get_z(wing.sections(index1).stringers(index2,x_pos),1);
        wing.sections(index1).stringers(index2,str_area) = 1; 
    end
end

%% Continue to follow code that Toohey gave us. 