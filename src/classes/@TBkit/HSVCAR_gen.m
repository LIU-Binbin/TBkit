function HSVCAR = HSVCAR_gen(orb_list,mode,discrimination,center,orientation)
%HSVCAR_GEN Generate HSV color values for orbital visualization
%
%   Syntax:
%       HSVCAR = HSVCAR_gen(orb_list,mode,discrimination,center,orientation)
%
%   Description:
%       Computes Hue-Saturation-Value (HSV) color representations for orbitals
%       based on their positions and selected visualization mode.
%
%   Inputs:
%       orb_list      - N×3 array of orbital positions
%       mode          - Visualization mode:
%                       'hinge', 'surf', 'select', 'select-points', 'orient', 'slab'
%       discrimination- Sensitivity parameter (default=0.1)
%       center        - Center/reference point(s) (default=[0.5,0.5,0.5])
%       orientation   - Directional axis (1=x,2=y,3=z, default=3)
%
%   Output:
%       HSVCAR - N×3 array of HSV color values (hue in first column)
%
%   See also: orbone2hsv, COLORCAR_gen
%--------  nargin  --------
if nargin < 3
    discrimination = 0.1;
end
if nargin < 4
    center  = [0.5, 0.5,0.5];
end

if nargin < 5
    orientation = 3;
end
if nargin < 2
    mode = 'hinge';
end
[norb,~] = size(orb_list);
HSVCAR  = zeros(norb,1);
switch mode
    case 'hinge'
        for i = 1:norb
            orb_one = orb_list(i,:);
            [~,~,HSVCAR(i)] = TBkit.orbone2hsv(orb_one,discrimination,center,orientation);
        end
    case 'surf'
        for i = 1:norb
            orb_one = orb_list(i,:);
            [~,HSVCAR(i),~] = TBkit.orbone2hsv(orb_one,discrimination,center,orientation);
        end
    case 'select'
        nCenter = size(center,1);
        for i = 1:nCenter/2
            %((abs(x)+1i*discrimination-0.5)*(abs(y)+1i*discrimination-0.5))^-1;
            HSVCAR = HSVCAR + ...
                imag(((abs(orb_list(:,1) -center(i,1))+1i*discrimination)...
                .*(abs(orb_list(:,2) -center(i,2))+1i*discrimination)...
                .*(abs(orb_list(:,3) -center(i,3))+1i*discrimination))...
                .^-1);
        end
        for i = nCenter/2+1 :nCenter
            %((abs(x)+1i*discrimination-0.5)*(abs(y)+1i*discrimination-0.5))^-1;
            HSVCAR = HSVCAR - ...
                imag(((abs(orb_list(:,1) -center(i,1))+1i*discrimination)...
                .*(abs(orb_list(:,2) -center(i,2))+1i*discrimination)...
                .*(abs(orb_list(:,3) -center(i,3))+1i*discrimination))...
                .^-1);
        end
        HSVCAR = HSVCAR/nCenter;
    case 'select-points'
        nCenter = size(center,1);
        for i = 1:nCenter
            %((abs(x)+1i*discrimination-0.5)*(abs(y)+1i*discrimination-0.5))^-1;
            HSVCAR = HSVCAR + ...
                imag(((abs(orb_list(:,1) -center(i,1))+1i*discrimination)...
                .*(abs(orb_list(:,2) -center(i,2))+1i*discrimination)...
                .*(abs(orb_list(:,3) -center(i,3))+1i*discrimination))...
                .^-1);
        end
        HSVCAR = HSVCAR/nCenter;
    case 'orient'
        for i = 1:norb
            orb_one = orb_list(i,:);
            [HSVCAR(i),~,~] = TBkit.orbone2hsv(orb_one,discrimination,center,orientation);
        end
    case 'slab'
        G_slab = abs((abs(orb_list(:,orientation)-center(orientation))-0.5)) < discrimination;
        HSVCAR = ~(G_slab);
end
% 0-1
WAN_NUM = norb;
HSVCAR = (-normalize(HSVCAR,'range')+1)/2;
HSVCAR(:,2:3) = ones(WAN_NUM,2);
end