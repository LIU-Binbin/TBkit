function RotationMat = nThetad2Rotation(thetad,n)
if nargin < 2
n = [0,0,1];
end
n_x = n(1);
n_y = n(2);
n_z = n(3);
RotationMat(1,1)=n_x^2*(1-cosd(thetad))+cosd(thetad);
RotationMat(1,2)=n_x*n_y*(1-cosd(thetad))-n_z*sind(thetad);
RotationMat(1,3)=n_x*n_z*(1-cosd(thetad))+n_y*sind(thetad);
RotationMat(2,1)=n_y*n_x*(1-cosd(thetad))+n_z*sind(thetad);
RotationMat(2,2)=n_y^2*(1-cosd(thetad))+cosd(thetad);
RotationMat(2,3)=n_y*n_z*(1-cosd(thetad))-n_x*sind(thetad);
RotationMat(3,1)=n_z*n_x*(1-cosd(thetad))-n_y*sind(thetad);
RotationMat(3,2)=n_z*n_y*(1-cosd(thetad))+n_x*sind(thetad);
RotationMat(3,3)=n_z^2*(1-cosd(thetad))+cosd(thetad);
end
