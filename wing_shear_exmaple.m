%wing shear flow
clear all;
close all;

% initial load:
Vx = 1; Vz = 0;

%define a few 
numTopStringers = 5;
numBottomStringers = 5;
numNoseTopStringers = 3;
numNoseBottomStringers = 3;

t_upper = 0.02/12;
t_lower = 0.02/12;
t_upper_front = 0.02/12;
t_lower_front = 0.02/12;

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

sparCaps(1).area = 1;
sparCaps(2).area = 1;
sparCaps(3).area = 1;
sparCaps(4).area = 1;

upperStringerGap = (sparCaps(3).posX - sparCaps(1).posX)/(numTopStringers + 1);
lowerStringerGap = (sparCaps(3).posX - sparCaps(1).posX)/(numBottomStringers + 1);
upperNoseStringerGap = (sparCaps(1).posX - 0)/(numNoseTopStringers + 1);
lowerNoseStringerGap = (sparCaps(1).posX - 0)/(numNoseBottomStringers + 1);


%assume stringers spaced evenly along X axis betwen Spars
%top Stringers
for i=1:numTopStringers
    topStringers(i).posX = sparCaps(1).posX + upperStringerGap*i;
    topStringers(i).posZ = get_z(topStringers(i).posX,1);
    topStringers(i).area = 1;
end

%bottom Stringers
for i=1:numBottomStringers
    bottomStringers(i).posX = sparCaps(1).posX + lowerStringerGap*i;
    bottomStringers(i).posZ = get_z(bottomStringers(i).posX,0);
    bottomStringers(i).area = 1;

end

%nose top Stringers
for i=1:numNoseTopStringers
    noseTopStringers(i).posX = upperNoseStringerGap*i;
    noseTopStringers(i).posZ = get_z(noseTopStringers(i).posX,1);
    noseTopStringers(i).area = 1;
end

%nose bottom Stringers
for i=1:numNoseBottomStringers
    noseBottomStringers(i).posX = lowerNoseStringerGap*i;
    noseBottomStringers(i).posZ = get_z(noseBottomStringers(i).posX,0);
    noseBottomStringers(i).area = 1;
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
for i=1:numBottomStringers %top stringers
    Ix = Ix + bottomStringers(i).area*(bottomStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + bottomStringers(i).area*(bottomStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + bottomStringers(i).area*(bottomStringers(i).posX-centroid.posX)*(bottomStringers(i).posZ-centroid.posZ);
end
for i=1:numNoseTopStringers %top stringers
    Ix = Ix + noseTopStringers(i).area*(noseTopStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + noseTopStringers(i).area*(noseTopStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + noseTopStringers(i).area*(noseTopStringers(i).posX-centroid.posX)*(noseTopStringers(i).posZ-centroid.posZ);
end
for i=1:numNoseBottomStringers %top stringers
    Ix = Ix + noseBottomStringers(i).area*(noseBottomStringers(i).posZ-centroid.posZ)^2;
    Iz = Iz + noseBottomStringers(i).area*(noseBottomStringers(i).posX-centroid.posX)^2;
    Ixz = Ixz + noseBottomStringers(i).area*(noseBottomStringers(i).posX-centroid.posX)*(noseBottomStringers(i).posZ-centroid.posZ)^2;
end

%define webs


%upper webs
for i=1:(numTopStringers+1)
    web(i).xStart = sparCaps(1).posX + upperStringerGap*(i-1);
    web(i).xEnd = sparCaps(1).posX + upperStringerGap*(i);
    web(i).thickness = t_upper;
    web(i).zStart = get_z(web(i).xStart,1);
    web(i).zEnd = get_z(web(i).xEnd,1);
    if i==1
        web(i).dp_area = sparCaps(1).area;
        web(i).dP = 0;
        web(i).qPrime = 0;
    else
        web(i).dp_area = topStringers(i-1).area;
        web(i).dP = get_dp(web(i).xStart-centroid.posX,web(i).zStart-centroid.posZ, ...
        Vx,Vz,Ix,Iz,Ixz,web(i).dp_area);
        web(i).qPrime = web(i).qPrime + web(i).dP;
    end
    
    tempInt = get_int(web(i).xStart,web(i).xEnd,1);  %integral of airfoil function
    triangle1 = (web(i).xStart - sparCaps(1).posX)*web(i).zStart/2;
    triangle2 = (web(i).xEnd - sparCaps(1).posX)*web(i).zEnd/2;
    web(i).Area = tempInt + triangle1 - triangle2;
    %web(i).dS = ...        Example:  get_distance(web(i).xStart,web(i).xEnd,1)
    %web(i).dS_over_t = ....
    %web(i).q_dS_over_t = ...
    
    %web(i).radCurv = ...   Example:  get_curve(web(i).xStart,web(i).xEnd,1)
end

% define more webs for spars and lower surfaces, and for both cells in
% order to solve for shear flow





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







