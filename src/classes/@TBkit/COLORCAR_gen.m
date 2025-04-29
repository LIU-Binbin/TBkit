function [COLORCAR,WEIGHTCAR] = COLORCAR_gen(WAVECAR,HSVCAR,signlist)
%COLORCAR_GEN Generate color and weight arrays for band visualization
%
%   Syntax:
%       [COLORCAR,WEIGHTCAR] = COLORCAR_gen(WAVECAR,HSVCAR,signlist)
%
%   Description:
%       Creates color and weight arrays for visualizing band structure
%       with orbital composition information.
%
%   Inputs:
%       WAVECAR  - Wavefunction array [orbitals×bands×k-points]
%       HSVCAR   - HSV color values for orbitals
%       signlist - Sign structure for orbitals
%
%   Outputs:
%       COLORCAR - RGB color structure array
%       WEIGHTCAR - Weight array for visualization
if nargin < 3
    signlist = 1;
end
if all(signlist==1)
    SIGN_mode = false;
else
    SIGN_mode = true;
end
[~,norb,kn] = size(WAVECAR);
temp.rgb = [0 0 0];
COLORCAR = repmat(temp,norb,kn);
WEIGHTCAR = zeros(norb,kn);
if SIGN_mode
    for ki = 1:kn
        WAVECAR_one = WAVECAR(:,:,ki);
        for orbi = 1:norb
            WAVEFUNC = WAVECAR_one(:,orbi);
            [WEIGHTCAR(orbi,ki),~]  = TBkit.COLOR_one_gen(WAVEFUNC,HSVCAR,signlist);
        end
    end
else
    for ki = 1:kn
        WAVECAR_one = WAVECAR(:,:,ki);
        for orbi = 1:norb
            WAVEFUNC = WAVECAR_one(:,orbi);
            [rgb,SIGN_one]  = TBkit.COLOR_one_gen(WAVEFUNC,HSVCAR);
            COLORCAR(orbi,ki).rgb = rgb;
            hsv_temp = rgb2hsv(rgb);
            WEIGHTCAR(orbi,ki) = SIGN_one*(hsv_temp(1)*100+hsv_temp(3)*10+hsv_temp(2)*1);
        end
    end
end
end