%%%%Arclength Calculations 
% Write Function for wing surfaces 
h=.001;
xChord = 0:h:1;

upperSurface = zeros(1,length(xChord));
lowerSurface = zeros(1,length(xChord));

for i=1:length(xChord)
    upperSurface(i) = get_z(xChord(i),1);
    lowerSurface(i) = get_z(xChord(i),0);
end

% Take first derivative
upperSurfacePrime = diff(upperSurface)/h;         
lowerSurfacePrime = diff(lowerSurface)/h;  

%setup integrand for arclength equation 
fupper = 1+upperSurfacePrime.^2;
flower = 1+lowerSurfacePrime.^2;

% Take integral 
Fupper=cumtrapz(fupper);
Flower=cumtrapz(flower);

%Use fact that integral from a to b of f(x) equals F(b)-F(a)
for index=1:numstingers 
    arclength = Fupper(StringerXpos(i+1)./(h))- Fupper(StringerXpos(i)./(h));
end
for index=1:numstingers 
    arclength = Flower(StringerXpos(i+1)./(h))- Flower(StringerXpos(i)./(h));
end
% % Finds arclength of between stringers fro a single section 
%%integrate upper surface


%%%%Calculating the area 
%integrate upper surface
intupper= cumtrapz(upperSurface);
%plug into area equation from slides
for index=1:numstingers 
    Area=uppersurface(stringerXpos(index))*(stringerXpos(index)-stringerXpos(index-1))/2-uppersurface(stringerXpos(index+1))*(tringerXpos(index+1)-stringerXpos(index-1))/2+(intupper(StringerXpos(i+1)./(h))- intupper(StringerXpos(i)./(h)));
end