%wing shear flow
clear all;
close all;
% 
% for wing_index=1:1000
%     
%     wings(wing_index).numTopStringers = round(rand*10)+10;
%     wings(wing_index).numBottomStringers = round(rand*10)+10;
%     wings(wing_index).frontSpar = 0.2;
%     wings(wing_index).backSpar = 0.75;
%     
%     
%     frontSpar = 0.2;
%     backSpar = 0.75;
% 
% 
%     sparCaps(1).posX = frontSpar;
%     sparCaps(2).posX = frontSpar;
%     sparCaps(3).posX = backSpar;
%     sparCaps(4).posX = backSpar;
%     
%     sparCaps(1).posZ = get_z(frontSpar,1);
%     sparCaps(2).posZ = get_z(frontSpar,0);
%     sparCaps(3).posZ = get_z(backSpar,1);
%     sparCaps(4).posZ = get_z(backSpar,0);
%     
%     wings(wing_index).sparCaps = sparCaps;
%     
%     %do all analysis
%     
%     %failure criteria
%     
% 
%     %is good?
%     
%     
%     
% end

%define a few 
numTopStringers = 10;
numBottomStringers = 2;
numNoseTopStringers = 2;
numNoseBottomStringers = 2;

frontSpar = 0.2;
backSpar = 0.75;

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

%Ix = ...
%Iz = ...
%Ixz= ...


%define webs


%upper webs
for i=1:length(numTopStringers)+1
    web(i).xStart = sparCaps(1).posX + upperStringerGap*(i-1);
    web(i).xEnd = sparCaps(1).posX + upperStringerGap*i;
    %web(i).thickness = ...
    %web(i).dp = ...  dP equation
    %web(i).qPrime = ...  summing dP
    %web(i).Area = ...      Example:  get_area(web(i).xStart,web(i).xEnd,1)
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







