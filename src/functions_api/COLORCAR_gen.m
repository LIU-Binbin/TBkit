function [COLORCAR, WEIGHTCAR] = COLORCAR_gen(WAVECAR, HSVCAR)
% COLORCAR_GEN Generates a color and weight matrix based on the wavefunctions (WAVECAR) and HSV color space (HSVCAR).
%
% Input:
%   WAVECAR  - A 3D matrix representing wavefunctions, with dimensions [n, norb, kn],
%              where n is the number of spatial points, norb is the number of orbitals, and kn is the number of k-points.
%   HSVCAR   - A matrix of HSV color representations for each orbital.
%
% Output:
%   COLORCAR - A structure array of RGB colors corresponding to each orbital and k-point.
%   WEIGHTCAR - A matrix of weights derived from the HSV components of the RGB color.

    % Get the size of the WAVECAR matrix
    [~, norb, kn] = size(WAVECAR);

    % Initialize the COLORCAR structure array and WEIGHTCAR matrix
    COLORCAR = struct('rgb', cell(norb, kn)); % Preallocate with empty cell arrays
    WEIGHTCAR = zeros(norb, kn); % Preallocate for weights

    % Loop over each k-point and orbital to compute the color and weight
    for ki = 1:kn
        WAVECAR_one = WAVECAR(:,:,ki); % Extract the wavefunction for this k-point

        for orbi = 1:norb
            WAVEFUNC = WAVECAR_one(:,orbi); % Extract the wavefunction for this orbital

            % Generate the RGB color for the wavefunction
            rgb = COLOR_one_gen(WAVEFUNC, HSVCAR);

            % Store the RGB color in the structure
            COLORCAR(orbi, ki).rgb = rgb;

            % Convert the RGB color to HSV and calculate the weight
            hsv_temp = rgb2hsv(rgb); % Convert RGB to HSV
            WEIGHTCAR(orbi, ki) = hsv_temp(1) * 100 + hsv_temp(3) * 10 + hsv_temp(2) * 1; % Calculate weight
        end
    end
end

% function [COLORCAR,WEIGHTCAR] = COLORCAR_gen(WAVECAR,HSVCAR)
%     [~,norb,kn] = size(WAVECAR);
%     temp.rgb = [0 0 0];
%     COLORCAR = repmat(temp,norb,kn);
%     WEIGHTCAR = zeros(norb,kn);
%     for ki = 1:kn
%         WAVECAR_one = WAVECAR(:,:,ki);
%         for orbi = 1:norb
%             WAVEFUNC = WAVECAR_one(:,orbi);
%             rgb = COLOR_one_gen(WAVEFUNC,HSVCAR);
%             COLORCAR(orbi,ki).rgb = rgb;
%             hsv_temp = rgb2hsv(rgb);
%             WEIGHTCAR(orbi,ki) = hsv_temp(1)*100+hsv_temp(3)*10+hsv_temp(2)*1;
%         end
%     end
% end