function [n,thetad]= Rotation2nThetad(R,Rm)
arguments
R double {mustBeReal,Oper.mustBeSize(R,[3,3])}
Rm = [1 0 0;0 1 0;0 0 1];
end
[n,theta]= Oper.Rotation2nTheta(R,Rm);
thetad = str2sym(Oper.name_angle(theta, false));
end
