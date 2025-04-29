function [Hue,surf_level,hing_level] = orbone2hsv(orb_one,discrimination,center,orientation)
%ORBONE2HSV Convert single orbital position to HSV color components
%
%   Syntax:
%       [Hue,surf_level,hing_level] = orbone2hsv(orb_one,discrimination,center,orientation)
%
%   Description:
%       Maps orbital position to color components using complex analysis
%       techniques. Used for visualizing orbital spatial distributions.
%
%   Inputs:
%       orb_one      - Orbital position [x,y,z]
%       discrimination - Sensitivity parameter (default=0.1)
%       center       - Reference center point (default=[0.5,0.5,0.5])
%       orientation  - Directional axis (1=x,2=y,3=z, default=3)
%
%   Outputs:
%       Hue        - Hue component (0-1)
%       surf_level - Surface detection level
%       hing_level - Hinge detection level
%
%   See also: HSVCAR_gen
%--------  nargin  --------
if nargin < 2
    discrimination = 0.1;
end
if nargin <3
    center  = [0.5, 0.5,0.5];
end
if nargin <4
    orientation = 3;
end

orb_init = orb_one - center;
%--------  init  --------
switch orientation
    case {1,-1}
        x = orb_init(2);
        y = orb_init(3);
    case {2,-2}
        x = orb_init(3);
        y = orb_init(1);
    case {3,-3}
        x = orb_init(1);
        y = orb_init(2);
end
r = norm([x y]);
z = x+1i*y;
Hue = (angle(z)+pi)/(2*pi);
G_hinge = ((abs(x)+1i*discrimination-0.5)*(abs(y)+1i*discrimination-0.5))^-1;
%G_hinge = (r+1i*discrimination-0.5*1.414)^-1;
G_surf = (r+1i*discrimination-0.5)^-1;
if orientation > 0
    hing_level = -imag(G_hinge);
    surf_level = -imag(G_surf);
else
    hing_level = -imag(G_hinge)*sign(x*y);
    surf_level = -imag(G_surf)*sign(x*y);
end

end