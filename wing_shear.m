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
num_spars       = 2; 

web             = struct('areas',0,'thickness',0.0017,'dp_area',0,'dP',0,'qPrime',0,'ds',0,'dS_over_t',0,'q_dS_over_t',0,'two_A_qprime',0,'qp_dx',0,'qp_dz',0);
webs            = web; 
for index = 2:sum(num_stringers)+num_sections-1+num_spars
    webs(index) = web;
end

section1    = struct('num_str', 0,'start_pos', spar1.position(1), 'end_pos', spar2.position(1), 'x_length', 0, 'stringers', zeros(num_stringers(1),3));
section2    = struct('num_str', 0,'start_pos', spar1.position(1), 'end_pos', spar2.position(1), 'x_length', 0, 'stringers', zeros(num_stringers(2),3)); 
section3    = struct('num_str', 0,'start_pos', 0                , 'end_pos', spar1.position(1), 'x_length', 0, 'stringers', zeros(num_stringers(3),3)); 
section4    = struct('num_str', 0,'start_pos', 0                , 'end_pos', spar1.position(1), 'x_length', 0, 'stringers', zeros(num_stringers(4),3));

% Variables to access stringer arrays. Faster than structures. 
x_pos       = 1; 
z_pos       = 2; 
str_area    = 3; 

%stringer = [x_position, y_position, area]
stringer    = [0,0,0]; 

%web      = start_x_position, end_x_position, thickness
web         = [0,0,0]; 

wing        = struct('spars',[spar1,spar2,spar3,spar4], 'sections', [section1, section2, section3, section4],'webs',webs);


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


%% Centroid

centroid_x  = 0;
centroid_x_area_sum=0;
centroid_z=0;
centroid_z_area_sum=0;
for index1 = 1:num_sections
    for index2 = 1:num_stringers(index1)
        centroid_x=centroid_x+(wing.sections(index1).stringers(index2,x_pos)...
            *wing.sections(index1).stringers(index2,str_area));
        centroid_x_area_sum=centroid_x_area_sum+wing.sections(index1).stringers(index2,str_area);
        
        centroid_z=centroid_z+(wing.sections(index1).stringers(index2,z_pos)...
            *wing.sections(index1).stringers(index2,str_area));
        centroid_z_area_sum=centroid_z_area_sum+wing.sections(index1).stringers(index2,str_area);
    end
end
for index1 = 1:num_spars
    centroid_x = centroid_x + wing.spars(index1).position(x_pos)*wing.spars(index1).area;
    centroid_x_area_sum = centroid_x_area_sum + wing.spars(index1).area;
    
    centroid_z = centroid_z + wing.spars(index1).position(z_pos)*wing.spars(index1).area;
    centroid_z_area_sum = centroid_z_area_sum + wing.spars(index1).area;
end


centroid_x = centroid_x/centroid_x_area_sum;
centroid_z=centroid_z/centroid_z_area_sum;

%% Moment of Inertia
Ix=0;Iz = 0; Ixz = 0; 
for index1 = 1:num_sections
    for index2 = 1:num_stringers(index1)
        Ix=Ix+(wing.sections(index1).stringers(index2,str_area)*...
            (wing.sections(index1).stringers(index2,z_pos)-centroid_z)^2);
        Iz=Iz+(wing.sections(index1).stringers(index2,str_area)*...
            (wing.sections(index1).stringers(index2,x_pos)-centroid_x)^2);
        Ixz = Ixz + wing.sections(index1).stringers(index2,str_area) *...
            (wing.sections(index1).stringers(index2,z_pos)-centroid_z) *...
            (wing.sections(index1).stringers(index2,x_pos)-centroid_x);
    end
end
for index1 = 1:num_spars
    Ix = Ix + wing.spars(index1).area*(wing.spars(index1).position(z_pos)-centroid_z)^2; 
    Iz = Iz + wing.spars(index1).area*(wing.spars(index1).position(x_pos)-centroid_x)^2; 
    Ixz = Ixz + wing.spars(index1).area*(wing.spars(index1).position(x_pos)-centroid_x)*...
                  (wing.spars(index1).position(z_pos)-centroid_z); 
end 
%% Web Calculations 

%Sample Inputs
Vx              = 1; 
Vz              = 0; 
index3          = 1; 
x_plus          = 0; 
x_minus         = 0; 
x_i             = 0;
z_i             = 0; 

%% Nested loop to run through all webs. The loops are used for all stringer
% webs and then there are conditions that handle all of the spar webs, and 
% differences between cells. 

% to update:
% 1) Make two web structs for cell1 and cell2. Break up code in here for
% different cells.
% 2) Update iterations for taking care if appropriate ordering of webs. 
% 3) Proof read all code. 

for index1 = 1:num_sections % run through sections
    for index2 = 1:(num_stringers(index1)+1) % run through webs. 
        
        % define positions for stringers to calculate areas. 
        if index2 == 1
            x_plus  = wing.sections(index1).stringers(index2);
            x_minus = wing.sections(index1).start_pos;
            x_i     = wing.sections(index1).start_pos; 
            z_i     = get_z(x_i, mod(index1,2)); 
            z_plus  = get_z(x_plus, mod(index1,2));
            
            integral = get_int(x_i,x_plus,mod(index1,2));
        elseif index2 == (num_stringers(index1)+1)
            x_plus  = wing.sections(index1).end_pos;
            x_minus = wing.sections(index1).start_pos;
            x_i     = wing.sections(index1).stringers(index2-1); 
            z_i     = get_z(x_i, mod(index1,2)); 
            z_plus  = get_z(x_plus, mod(index1,2));
            integral = get_int(x_i,x_plus,mod(index1,2));
        else
            x_plus  = wing.sections(index1).stringers(index2);
            x_minus = wing.sections(index1).start_pos;
            x_i     = wing.sections(index1).stringers(index2-1); 
            z_i     = get_z(x_i, mod(index1,2)); 
            z_plus  = get_z(x_plus, mod(index1,2));
            integral = get_int(x_i,x_plus,mod(index1,2));
        end
        
        % Sum areas from triangles and integrals above. 
        wing.webs(index3).areas = abs(z_i*(x_i-x_minus)/2) - abs((z_plus*(x_plus-x_minus))/2) + integral(end);     
        
        % outline special cases for shear flow starting conditions, based
        % off of cell and section. 
        if index1 == 1 && index2 == 1
            wing.webs(index3).dp_area = wing.spars(index1).area; 
            wing.webs(index3).dP      = 0; 
            wing.webs(index3).qPrime  = 0; 
        elseif index1 == 2 && index2 == 1
            wing.webs(index3).dp_area = wing.spars(3).area; 
            %update this for changing shear flow connection points. 
            wing.webs(index3).dP      = get_dp(x_i - centroid_x,z_i-centroid_z,Vx,Vz,Ix,Iz,Ixz,wing.webs(index3).dp_area); 
            wing.webs(index3).qPrime  = 0; %% input value from rear spar here. 
        elseif index1 == 3 && index2 == 1
            wing.webs(index3).dp_area = wing.spars(4).area; 
            wing.webs(index3).dP      = 0;
            wing.webs(index3).qPrime  = 0;
        elseif index1 == 4 && index2 == 1
            wing.webs(index3).dp_area = 0; 
            % change this for altering shear flow connection points. 
            wing.webs(index3).dP      = get_dp(x_i - centroid_x,z_i-centroid_z,Vx,Vz,Ix,Iz,Ixz,wing.webs(index3).dp_area); 
            wing.webs(index3).qPrime  = 0; %% input value from lower nose here. 
        else
            wing.webs(index3).dp_area = wing.sections(index1).stringers(index2-1,str_area); 
            wing.webs(index3).dP      = get_dp(x_i - centroid_x, z_i-centroid_z,Vx,Vz,Ix,Iz,Ixz,wing.webs(index3).dp_area); 
            wing.webs(index3).qPrime  = wing.webs(index3-1).qPrime - wing.webs(index3).dP;   
        end
        
        wing.webs(index3).ds            = get_ds(x_i,x_plus,mod(index1,2)); 
        wing.webs(index3).dS_over_t     = wing.webs(index3).ds / wing.webs(index3).thickness;
        wing.webs(index3).q_dS_over_t   = wing.webs(index3).qPrime*wing.webs(index3).dS_over_t; 
        wing.webs(index3).two_A_qprime  = 2*wing.webs(index3).areas*wing.webs(index3).qPrime; 
        wing.webs(index3).qp_dx         = wing.webs(index3).qPrime*(x_plus - x_i); 
        wing.webs(index3).qp_dz         = wing.webs(index3).qPrime*(z_plus - z_i); 

        index3 = index3+1;

    end
end

%% Plots
% 
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
