function RotationMat = nTheta2Rotation(theta,n)
if nargin < 2
n = [0,0,1];
end
n = normalize(n);
thetad = rad2deg(theta);
RotationMat = Oper.nThetad2Rotation(thetad,n);
end
