%wing shear flow
clear all;
close all;


for testLoop =1:2

% initial load:
if testLoop == 2
    Vx = 1; Vz = 0;
else 
    Vx = 0; Vz = 1;    
end

%define a few 
numTopStringers = 5;
numBottomStringers = 5;
numNoseTopStringers = 3;
numNoseBottomStringers = 3;

t_upper = 0.02/12;
t_lower = 0.02/12;
t_upper_front = 0.02/12;
t_lower_front = 0.02/12;
t_frontSpar = 0.04/12;
t_rearSpar = 0.04/12;

frontSpar = 0.2;
backSpar = 0.7;

sparCaps(1).posX = frontSpar;
sparCaps(2).posX = frontSpar;
sparCaps(3).posX = backSpar;
sparCaps(4).posX = backSpar;

sparCaps(1).posZ = get_z(frontSpar,1);
sparCaps(2).posZ = get_z(frontSpar,0);
sparCaps(3).posZ = get_z(backSpar,1);
sparCaps(4).posZ = get_z(backSpar,0);

sparCaps(1).area = .1;
sparCaps(2).area = .1;
sparCaps(3).area = .1;
sparCaps(4).area = .1;

upperStringerGap = (sparCaps(3).posX - sparCaps(1).posX)/(numTopStringers + 1);
lowerStringerGap = (sparCaps(3).posX - sparCaps(1).posX)/(numBottomStringers + 1);
upperNoseStringerGap = (sparCaps(1).posX - 0)/(numNoseTopStringers + 1);
lowerNoseStringerGap = (sparCaps(1).posX - 0)/(numNoseBottomStringers + 1);


%assume stringers spaced evenly along X axis betwen Spars
%top Stringers
for i=1:numTopStringers
    topStringers(i).posX = sparCaps(1).posX + upperStringerGap*i;
    topStringers(i).posZ = get_z(topStringers(i).posX,1);
    topStringers(i).area = .1;
end

%bottom Stringers
for i=1:numBottomStringers
    bottomStringers(i).posX = sparCaps(4).posX - lowerStringerGap*i;
    bottomStringers(i).posZ = get_z(bottomStringers(i).posX,0);
    bottomStringers(i).area = .1;

end

%nose bottom Stringers
for i=1:numNoseBottomStringers
    noseBottomStringers(i).posX = sparCaps(2).posX - lowerNoseStringerGap*i;
    noseBottomStringers(i).posZ = get_z(noseBottomStringers(i).posX,0);
    noseBottomStringers(i).area = .1;
end

%nose top Stringers
for i=1:numNoseTopStringers
    noseTopStringers(i).posX = upperNoseStringerGap*i;
    noseTopStringers(i).posZ = get_z(noseTopStringers(i).posX,1);
    noseTopStringers(i).area = .1;
end


centroid.posX = sum([sparCaps.posX].*[sparCaps.area]) + ...
    sum([topStringers.posX].*[topStringers.area]) + ...
    sum([bottomStringers.posX].*[bottomStringers.area]) + ...
    sum([noseTopStringers.posX].*[noseTopStringers.area]) + ...
    sum([noseBottomStringers.posX].*[noseBottomStringers.area]);

centroid.posX = centroid.posX / ( sum([sparCaps.area]) + sum([topStringers.area]) + ...
    sum([bottomStringers.area]) + sum([noseTopStringers.area]) + sum([noseBottomStringers.area]));

centroid.posZ = sum([sparCaps.posZ].*[sparCaps.area]) + ...
    sum([topStringers.posZ].*[topStringers.area]) + ...
    sum([bottomStringers.posZ].*[bottomStringers.area]) + ...
    sum([noseTopStringers.posZ].*[noseTopStringers.area]) + ...
    sum([noseBottomStringers.posZ].*[noseBottomStringers.area]);
centroid.posZ = centroid.posZ / ( sum([sparCaps.area]) + sum([topStringers.area]) + ...
    sum([bottomStringers.area]) + sum([noseTopStringers.area]) + sum([noseBottomStringers.area]));

%summing contributions for inertia terms
Ix = 0; Iz = 0; Ixz = 0;

for i=1:4 %spar caps
    Ix = Ix + sparCaps(i).area*(sparCaps(i).posZ-centroid.posZ)^2;
    Iz = Iz + sparCaps(i).area*(sparCaps(i).posX-centroid.posX)^2;
    Ixz = Ixz + sparCaps(i).area*(sparCaps(i).posX-centroid.posX)*(sparCaps(i).posZ-centroid.posZ);
end


for i=1:numTopStringers %top stringers
    Ix = Ix + topStringers(i).area*(topStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + topStringers(i).area*(topStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + topStringers(i).area*(topStringers(i).posX-centroid.posX)*(topStringers(i).posZ-centroid.posZ);
end
for i=1:numBottomStringers %bottom stringers
    Ix = Ix + bottomStringers(i).area*(bottomStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + bottomStringers(i).area*(bottomStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + bottomStringers(i).area*(bottomStringers(i).posX-centroid.posX)*(bottomStringers(i).posZ-centroid.posZ);
end
for i=1:numNoseTopStringers %nose top stringers
    Ix = Ix + noseTopStringers(i).area*(noseTopStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + noseTopStringers(i).area*(noseTopStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + noseTopStringers(i).area*(noseTopStringers(i).posX-centroid.posX)*(noseTopStringers(i).posZ-centroid.posZ);
end
for i=1:numNoseBottomStringers %nose bottom stringers
    Ix = Ix + noseBottomStringers(i).area*(noseBottomStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + noseBottomStringers(i).area*(noseBottomStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + noseBottomStringers(i).area*(noseBottomStringers(i).posX-centroid.posX)*(noseBottomStringers(i).posZ-centroid.posZ);
end

%Ixz = -Ixz;

%define webs

%% web cell 1

%upper webs
numStringers = numTopStringers;
stringerGap = upperStringerGap;
webThickness = t_upper;
tempStringers = topStringers;

for i=1:(numStringers+1)
    web(i).xStart = sparCaps(1).posX + stringerGap*(i-1);
    web(i).xEnd = sparCaps(1).posX + stringerGap*(i);
    web(i).thickness = webThickness;
    web(i).zStart = get_z(web(i).xStart,1);
    web(i).zEnd = get_z(web(i).xEnd,1);
    if i==1
        web(i).dp_area = sparCaps(1).area;
        web(i).dP = 0;
        web(i).qPrime = 0;
    else
        web(i).dp_area = tempStringers(i-1).area;
        web(i).dP = get_dp(web(i).xStart-centroid.posX,web(i).zStart-centroid.posZ, ...
        Vx,Vz,Ix,Iz,Ixz,web(i).dp_area);
        web(i).qPrime = web(i-1).qPrime - web(i).dP;
    end
    tempInt = get_int(web(i).xStart,web(i).xEnd,1);  %integral of airfoil function
    triangle1 = abs( (web(i).xStart - sparCaps(1).posX)*web(i).zStart/2);
    triangle2 = abs((web(i).xEnd - sparCaps(1).posX)*web(i).zEnd/2);
    web(i).Area = tempInt + triangle1 - triangle2;
    web(i).ds = get_ds(web(i).xStart,web(i).xEnd,1);
    web(i).dS_over_t = web(i).ds / web(i).thickness;
    web(i).q_dS_over_t = web(i).qPrime * web(i).dS_over_t;
    web(i).two_A_qprime = 2*web(i).Area*web(i).qPrime;
    web(i).qp_dx = web(i).qPrime*(web(i).xEnd-web(i).xStart);
    web(i).qp_dz = web(i).qPrime*(web(i).zEnd-web(i).zStart);
    %web(i).radCurv = ...   Example:  get_curve(web(i).xStart,web(i).xEnd,1)
end
webTop = web;
web = [];

%rear spar
i=1;
web(i).xStart = sparCaps(3).posX;
web(i).xEnd = sparCaps(4).posX;
web(i).thickness = t_rearSpar;
web(i).zStart = sparCaps(3).posZ;
web(i).zEnd = sparCaps(4).posZ;
web(i).dp_area = sparCaps(3).area;
web(i).dP = get_dp(web(i).xStart-centroid.posX,web(i).zStart-centroid.posZ, ...
        Vx,Vz,Ix,Iz,Ixz,web(i).dp_area);
web(i).qPrime = webTop(numTopStringers+1).qPrime - web(i).dP;
web(i).Area = (sparCaps(3).posX-sparCaps(1).posX)*sparCaps(3).posZ/2 + ...
    abs((sparCaps(3).posX-sparCaps(1).posX)*sparCaps(4).posZ/2);
web(i).ds = abs(sparCaps(3).posZ - sparCaps(4).posZ);
web(i).dS_over_t = web(i).ds / web(i).thickness;
web(i).q_dS_over_t = web(i).qPrime * web(i).dS_over_t;
web(i).two_A_qprime = 2*web(i).Area*web(i).qPrime;
    web(i).qp_dx = web(i).qPrime*(web(i).xEnd-web(i).xStart);
    web(i).qp_dz = web(i).qPrime*(web(i).zEnd-web(i).zStart);

webRearSpar = web;
web = [];


%lower webs
numStringers = numBottomStringers;
stringerGap = lowerStringerGap;
webThickness = t_lower;
tempStringers = bottomStringers;

for i=1:(numStringers+1)
    web(i).xStart = sparCaps(4).posX - stringerGap*(i-1);
    web(i).xEnd = sparCaps(4).posX - stringerGap*(i);
    web(i).thickness = webThickness;
    web(i).zStart = get_z(web(i).xStart,0);
    web(i).zEnd = get_z(web(i).xEnd,0);
    if i==1
        web(i).dp_area = sparCaps(4).area;
        web(i).dP = get_dp(web(i).xStart-centroid.posX,web(i).zStart-centroid.posZ, ...
        Vx,Vz,Ix,Iz,Ixz,web(i).dp_area);
        web(i).qPrime = webRearSpar.qPrime - web(i).dP;
    else
        web(i).dp_area = tempStringers(i-1).area;
        web(i).dP = get_dp(web(i).xStart-centroid.posX,web(i).zStart-centroid.posZ, ...
        Vx,Vz,Ix,Iz,Ixz,web(i).dp_area);
        web(i).qPrime = web(i-1).qPrime - web(i).dP;
    end
    
    tempInt = get_int(web(i).xEnd,web(i).xStart,0);  %integral of airfoil function
    triangle2 = abs((web(i).xStart - sparCaps(1).posX)*web(i).zStart/2);
    triangle1 = abs((web(i).xEnd - sparCaps(1).posX)*web(i).zEnd/2);
    web(i).Area = tempInt + triangle1 - triangle2;
    web(i).ds = get_ds(web(i).xStart,web(i).xEnd,0);
    web(i).dS_over_t = web(i).ds / web(i).thickness;
    web(i).q_dS_over_t = web(i).qPrime * web(i).dS_over_t;
    web(i).two_A_qprime = 2*web(i).Area*web(i).qPrime;
    web(i).qp_dx = web(i).qPrime*(web(i).xEnd-web(i).xStart);
    web(i).qp_dz = web(i).qPrime*(web(i).zEnd-web(i).zStart);

    %web(i).radCurv = ...   Example:  get_curve(web(i).xStart,web(i).xEnd,1)
end
webBottom = web;
web = [];

%front Spar
i=1;
web(i).xStart = sparCaps(2).posX;
web(i).xEnd = sparCaps(1).posX;
web(i).thickness = t_frontSpar;
web(i).zStart = sparCaps(2).posZ;
web(i).zEnd = sparCaps(1).posZ;
web(i).dp_area = sparCaps(2).area;
web(i).dP = get_dp(web(i).xStart-centroid.posX,web(i).zStart-centroid.posZ, ...
        Vx,Vz,Ix,Iz,Ixz,web(i).dp_area);
web(i).qPrime = webBottom(numBottomStringers+1).qPrime - web(i).dP;
web(i).Area = 0;
web(i).ds = abs(sparCaps(2).posZ - sparCaps(1).posZ);
web(i).dS_over_t = web(i).ds / web(i).thickness;
web(i).q_dS_over_t = web(i).qPrime * web(i).dS_over_t;
web(i).two_A_qprime = 2*web(i).Area*web(i).qPrime;
    web(i).qp_dx = web(i).qPrime*(web(i).xEnd-web(i).xStart);
    web(i).qp_dz = web(i).qPrime*(web(i).zEnd-web(i).zStart);

webFrontSpar = web;
web = [];




%% web cell 2

%lower nose webs
numStringers = numNoseBottomStringers;
stringerGap = lowerNoseStringerGap;
webThickness = t_lower_front;
tempStringers = noseBottomStringers;

for i=1:(numStringers+1)
    web(i).xStart = sparCaps(2).posX - stringerGap*(i-1);
    web(i).xEnd = sparCaps(2).posX - stringerGap*(i);
    web(i).thickness = webThickness;
    web(i).zStart = get_z(web(i).xStart,0);
    web(i).zEnd = get_z(web(i).xEnd,0);
    if i==1
        web(i).dp_area = sparCaps(2).area;
        web(i).dP = 0;
        web(i).qPrime = 0;
    else
        web(i).dp_area = tempStringers(i-1).area;
        web(i).dP = get_dp(web(i).xStart-centroid.posX,web(i).zStart-centroid.posZ, ...
        Vx,Vz,Ix,Iz,Ixz,web(i).dp_area);
        web(i).qPrime = web(i-1).qPrime - web(i).dP;
    end
    tempInt = get_int(web(i).xEnd,web(i).xStart,0);  %integral of airfoil function
    triangle1 = abs((web(i).xStart - sparCaps(2).posX)*web(i).zStart/2);
    triangle2 = abs((web(i).xEnd - sparCaps(2).posX)*web(i).zEnd/2);
    web(i).Area = tempInt + triangle1 - triangle2;
    web(i).ds = get_ds(web(i).xStart,web(i).xEnd,0);
    web(i).dS_over_t = web(i).ds / web(i).thickness;
    web(i).q_dS_over_t = web(i).qPrime * web(i).dS_over_t;
    web(i).two_A_qprime = 2*web(i).Area*web(i).qPrime;
    web(i).qp_dx = web(i).qPrime*(web(i).xEnd-web(i).xStart);
    web(i).qp_dz = web(i).qPrime*(web(i).zEnd-web(i).zStart);

    %web(i).radCurv = ...   Example:  get_curve(web(i).xStart,web(i).xEnd,1)
end
webLowerNose = web;
web = [];

%upper nose webs
numStringers = numNoseTopStringers;
stringerGap = upperNoseStringerGap;
webThickness = t_upper_front;
tempStringers = noseTopStringers;

for i=1:(numStringers+1)
    web(i).xStart = stringerGap*(i-1);
    web(i).xEnd = stringerGap*(i);
    web(i).thickness = webThickness;
    web(i).zStart = get_z(web(i).xStart,1);
    web(i).zEnd = get_z(web(i).xEnd,1);
    if i==1
        web(i).dp_area = 0;
        web(i).dP = 0;
        web(i).qPrime = webLowerNose(numNoseBottomStringers+1).qPrime - web(i).dP;
    else
        web(i).dp_area = tempStringers(i-1).area;
        web(i).dP = get_dp(web(i).xStart-centroid.posX,web(i).zStart-centroid.posZ, ...
        Vx,Vz,Ix,Iz,Ixz,web(i).dp_area);
        web(i).qPrime = web(i-1).qPrime - web(i).dP;
    end
    tempInt = get_int(web(i).xStart,web(i).xEnd,1);  %integral of airfoil function
    triangle2 = abs((web(i).xStart - sparCaps(2).posX)*web(i).zStart/2);
    triangle1 = abs((web(i).xEnd - sparCaps(2).posX)*web(i).zEnd/2);
    web(i).Area = tempInt + triangle1 - triangle2;
    web(i).ds = get_ds(web(i).xStart,web(i).xEnd,1);
    web(i).dS_over_t = web(i).ds / web(i).thickness;
    web(i).q_dS_over_t = web(i).qPrime * web(i).dS_over_t;
    web(i).two_A_qprime = 2*web(i).Area*web(i).qPrime;
    web(i).qp_dx = web(i).qPrime*(web(i).xEnd-web(i).xStart);
    web(i).qp_dz = web(i).qPrime*(web(i).zEnd-web(i).zStart);

    %web(i).radCurv = ...   Example:  get_curve(web(i).xStart,web(i).xEnd,1)
end
webUpperNose = web;
web = [];


%front Spar
i=1;
web(i).xStart = sparCaps(1).posX;
web(i).xEnd = sparCaps(2).posX;
web(i).thickness = t_frontSpar;
web(i).zStart = sparCaps(1).posZ;
web(i).zEnd = sparCaps(2).posZ;
web(i).dp_area = sparCaps(1).area;
web(i).dP = get_dp(web(i).xStart-centroid.posX,web(i).zStart-centroid.posZ, ...
        Vx,Vz,Ix,Iz,Ixz,web(i).dp_area);
web(i).qPrime = webUpperNose(numNoseTopStringers+1).qPrime - web(i).dP;
web(i).Area = 0;
web(i).ds = abs(sparCaps(1).posZ - sparCaps(2).posZ);
web(i).dS_over_t = web(i).ds / web(i).thickness;
web(i).q_dS_over_t = web(i).qPrime * web(i).dS_over_t;
web(i).two_A_qprime = 2*web(i).Area*web(i).qPrime;
    web(i).qp_dx = web(i).qPrime*(web(i).xEnd-web(i).xStart);
    web(i).qp_dz = web(i).qPrime*(web(i).zEnd-web(i).zStart);

webFrontSparCell2 = web;
web = [];


%check that dp*dx sums up to Vx

Fx = sum([webTop.qp_dx])+webRearSpar.qp_dx+ sum([webBottom.qp_dx])+webFrontSpar.qp_dx;  %cell 1
Fx = Fx + sum([webLowerNose.qp_dx])+ sum([webUpperNose.qp_dx]);  %cell 2
%Fx
Fz = sum([webTop.qp_dz])+webRearSpar.qp_dz+ sum([webBottom.qp_dz])+webFrontSpar.qp_dz;  %cell 1
Fz = Fz + sum([webLowerNose.qp_dz])+ sum([webUpperNose.qp_dz]);  %cell 2
%Fz
%%

% sum up the ds/t and  q*ds/t to solve 2 equations, 2 unknowns

% [A]*[q1s q2s] = B

A11 = sum([webTop.dS_over_t])+webRearSpar.dS_over_t+ sum([webBottom.dS_over_t])+webFrontSpar.dS_over_t;
A22 = sum([webLowerNose.dS_over_t])+ sum([webUpperNose.dS_over_t])+webFrontSparCell2.dS_over_t;
A12 = -webFrontSpar.dS_over_t;
A21 = -webFrontSparCell2.dS_over_t;

B1 = sum([webTop.q_dS_over_t])+webRearSpar.q_dS_over_t+ sum([webBottom.q_dS_over_t])+webFrontSpar.q_dS_over_t;
B2 = sum([webLowerNose.q_dS_over_t])+ sum([webUpperNose.q_dS_over_t])+webFrontSparCell2.q_dS_over_t;

Amat = [A11 A12; A21 A22];
Bmat = -[B1;B2];

qs = inv(Amat)*Bmat;



sum_2_a_q = sum([webTop.two_A_qprime])+webRearSpar.two_A_qprime+ sum([webBottom.two_A_qprime]);  %cell 1 qprimes
sum_2_a_q = sum_2_a_q + sum([webLowerNose.two_A_qprime])+ sum([webUpperNose.two_A_qprime]);   %cell 2 qprimes
sum_2_a_q = sum_2_a_q +  qs(1)*(sum([webTop.Area])+webRearSpar.Area+ sum([webBottom.Area]));
sum_2_a_q = sum_2_a_q +  qs(2)*(sum([webLowerNose.Area])+ sum([webUpperNose.Area]));

if Vx == 0
   sc.posX =  sum_2_a_q / Vz + frontSpar;
elseif Vz == 0 
   sc.posZ =  sum_2_a_q / Vx;
end


% now consider the torque representing shifting the load from the quarter
% chord to the SC  (need to check signs on these moments)
 
if Vx == 0
     torque = Vz*(sc.posX - 0.25);
elseif Vz == 0
    torque = Vx*sc.posZ;
end

Area1 = sum([webTop.Area]) + webRearSpar.Area + sum([webBottom.Area]);
%check area
Area1_check = get_int(frontSpar,backSpar,1) + get_int(frontSpar,backSpar,0);

Area2 = sum([webLowerNose.Area]) + sum([webUpperNose.Area]);
Area2_check = get_int(0,frontSpar,1) + get_int(0,frontSpar,0);


%for twist equation  (see excel spreadsheet example)

q1t_over_q2t = (A22/Area2 + webFrontSpar.dS_over_t/Area1)/(A11/Area1 + webFrontSpar.dS_over_t/Area2);

q2t = torque/(2*Area1*q1t_over_q2t + 2*Area2);
q1t = q2t*q1t_over_q2t;

qt = [q1t;q2t];

%--- insert force balance to check total shear flows --- 

% --- -- 


end

sc


%plotting airfoil cross-section

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

plot([sparCaps(1).posX sparCaps(2).posX],[sparCaps(1).posZ sparCaps(2).posZ],'-')
plot([sparCaps(3).posX sparCaps(4).posX],[sparCaps(3).posZ sparCaps(4).posZ],'-')
plot([sparCaps.posX],[sparCaps.posZ],'o')
plot([topStringers.posX],[topStringers.posZ],'or')
plot([bottomStringers.posX],[bottomStringers.posZ],'or')
plot([noseTopStringers.posX],[noseTopStringers.posZ],'or')
plot([noseBottomStringers.posX],[noseBottomStringers.posZ],'or')
plot(centroid.posX,centroid.posZ,'rx')
plot(sc.posX,sc.posZ,'gx')
