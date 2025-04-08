function COLOR_one = COLOR_one_gen(WF, HSVCAR)
% COLOR_ONE_GEN Generates a single color by weighting HSV values with the magnitude
% of the wavefunction (WF). The resulting color is normalized to the [0, 1] range.
%
% Input:
%   WF      - Wavefunction or weights, a vector of complex numbers.
%   HSVCAR  - A matrix where each row is an HSV color representation [H, S, V].
%
% Output:
%   COLOR_one - A 1x3 vector representing the resulting RGB color.

    % Calculate the magnitude squared of the wavefunction (WF)
    n = abs(WF).^2;

    % Convert the HSV colors to RGB
    RGBCAR = hsv2rgb(HSVCAR);

    % Calculate the weighted sum of RGB values based on WF magnitude
    COLOR_one = sum(RGBCAR .* n, 1);

    % Normalize the resulting color to the [0, 1] range
    COLOR_one = normalize(COLOR_one, 'range');
end


% function COLOR_one = COLOR_one_gen(WF,HSVCAR)
%     n = abs(WF.*conj(WF));
%     %[WAN_NUM,~] = size(HSVCAR);
%     %HSVCAR(:,2:3) = ones(WAN_NUM,2);
%     RGBCAR =  hsv2rgb(HSVCAR);
%     COLOR_one = sum(RGBCAR.*n,1);
%     COLOR_one = normalize(COLOR_one,'range');
% end