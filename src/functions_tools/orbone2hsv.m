function [Hue, surf_level, hing_level] = orbone2hsv(orb_one, discrimination, center, orientation)
    % orbone2hsv converts orb_one (coordinates) to HSV values.
    % Inputs:
    %   orb_one      - 1x3 vector of coordinates [x, y, z]
    %   discrimination - value to avoid zero division, default is 0.1
    %   center       - center of transformation, default is [0.5, 0.5, 0.5]
    %   orientation  - orientation for coordinate swapping, default is 3
    % Outputs:
    %   Hue          - Hue component of the HSV
    %   surf_level   - Surface level (related to distance in the plane)
    %   hing_level   - Hinge level (related to depth in the plane)

    % Set default values for optional inputs if not provided
    arguments
        orb_one 
        discrimination = 0.1;   % Default value for discrimination
        center = [0.5, 0.5, 0.5]; % Default center 
        orientation = 3;  % Default orientation
    end

    % Subtract the center from orb_one to get orb_init
    orb_init = orb_one - center;

    % Switch case for orientation and assigning x, y based on the selected orientation
    switch orientation
        case 1
            x = orb_init(2);  % Use y as x
            y = orb_init(3);  % Use z as y
        case 2
            x = orb_init(3);  % Use z as x
            y = orb_init(1);  % Use x as y
        case 3
            x = orb_init(1);  % Use x as x
            y = orb_init(2);  % Use y as y
    end

    % Calculate the radial distance (r) and complex number z
    r = norm([x, y]);
    z = x + 1i * y;  % Complex number representation for angle calculation

    % Calculate the Hue (angle of the complex number in the plane)
    Hue = (angle(z) + pi) / (2 * pi);  % Normalize angle to [0, 1]

    % Calculate hinge and surface levels using complex arithmetic
    G_hinge = ((abs(x) + 1i * discrimination - 0.5) * (abs(y) + 1i * discrimination - 0.5))^-1;
    G_surf = (r + 1i * discrimination - 0.5)^-1;

    % Extract the imaginary parts for hinge and surface levels
    hing_level = -imag(G_hinge);
    surf_level = -imag(G_surf);
end
