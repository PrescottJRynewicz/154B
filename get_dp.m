function z = get_dp(xDist,zDist,Vx,Vz,Ix,Iz,Ixz,A)

denom = (Ix*Iz-Ixz^2);
z = A*xDist*(Ix*Vx+Ixz*Vz)/denom - A*zDist*(Iz*Vz+Ixz*Vx)/denom;