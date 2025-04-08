function HSVCAR = HSVCAR_gen(orb_list, mode, discrimination, center, orientation)
% HSVCAR_gen - Generate the HSVCAR values based on orbital list and specified parameters.
%
% Usage:
%   HSVCAR = HSVCAR_gen(orb_list, mode, discrimination, center, orientation)
%
% Inputs:
%   orb_list      - A list of orbitals with their coordinates (each row is an orbital).
%   mode          - The mode of calculation, options are 'hinge', 'surf', and 'orient'.
%   discrimination - Discrimination parameter used for classification (default: 0.1).
%   center        - The reference center for orbital calculations (default: [0.5, 0.5, 0.5]).
%   orientation   - The orientation parameter affecting the calculation (default: 3).
%
% Outputs:
%   HSVCAR        - A vector of HSVCAR values computed for each orbital.

    % Input validation and default settings
    if nargin < 3
        discrimination = 0.1;  % Default discrimination value
    end
    if nargin < 4
        center = [0.5, 0.5, 0.5];  % Default center
    end
    if nargin < 5
        orientation = 3;  % Default orientation
    end
    if nargin < 2
        mode = 'hinge';  % Default mode is 'hinge'
    end

    % Get the number of orbitals
    [norb, ~] = size(orb_list);

    % Initialize HSVCAR with zeros
    HSVCAR = zeros(norb, 1);

    % Process based on the selected mode
    switch mode
        case 'hinge'
            % For 'hinge' mode, calculate the HSVCAR values
            for i = 1:norb
                orb_one = orb_list(i, :);
                [~, ~, HSVCAR(i)] = orbone2hsv(orb_one, discrimination, center, orientation);
            end
        case 'surf'
            % For 'surf' mode, calculate the HSVCAR values
            for i = 1:norb
                orb_one = orb_list(i, :);
                [~, HSVCAR(i), ~] = orbone2hsv(orb_one, discrimination, center, orientation);
            end
        case 'orient'
            % For 'orient' mode, calculate the HSVCAR values
            for i = 1:norb
                orb_one = orb_list(i, :);
                [HSVCAR(i), ~, ~] = orbone2hsv(orb_one, discrimination, center, orientation);
            end
        otherwise
            error('Invalid mode selected. Options are: "hinge", "surf", or "orient".');
    end

    % Normalize the HSVCAR values to the range [0, 1]
    HSVCAR = (-normalize(HSVCAR, 'range') + 1) / 2;

    % Add extra columns to HSVCAR matrix
    HSVCAR(:, 2:3) = ones(norb, 2);  % Set the second and third columns to ones
end
