function [COLOR_one,SIGN_one] = COLOR_one_gen(WF,HSVCAR,signlist)
%COLOR_ONE_GEN Generate color and sign for wavefunction visualization
%
%   Syntax:
%       [COLOR_one,SIGN_one] = COLOR_one_gen(WF,HSVCAR,signlist)
%
%   Description:
%       Computes color and sign for wavefunction plotting based on
%       orbital weights and sign structure.
%
%   Inputs:
%       WF       - Wavefunction vector
%       HSVCAR   - HSV color values for orbitals
%       signlist - Sign structure for orbitals
%
%   Outputs:
%       COLOR_one - Composite color value
%       SIGN_one  - Overall sign

if nargin <3
    SIGN_one = 1;
    signlist = 1;
else

end

n = abs(WF.*conj(WF));
%[WAN_NUM,~] = size(HSVCAR);
%HSVCAR(:,2:3) = ones(WAN_NUM,2);
if nargin == 3
    COLOR_one = sum(signlist.*HSVCAR(:,1).*n,1);
    %COLOR_one = normalize(COLOR_one,'range',[-1,1]);
else
    RGBCAR =  hsv2rgb(HSVCAR);
    COLOR_one = sum(RGBCAR.*n,1);
    COLOR_one = normalize(COLOR_one,'range');
end
if nargin == 3
    SIGN_one = sign(sum(signlist.*n,1));
end

if all(signlist==1)
    SIGN_one = 1;
    return;
end

end