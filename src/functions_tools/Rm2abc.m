function R_struct = Rm2abc(Rm)
    % Rm2abc converts a 3x3 rotation matrix Rm to its corresponding
    % 3 lattice parameters: a, b, c, and angles alpha, beta, gamma.
    % Input: Rm - 3x3 rotation matrix
    % Output: R_struct - structure containing lattice parameters
    
    % Compute the norms of the rows in one go
    a = norm(Rm(1,:));  % Length of vector a (row 1)
    b = norm(Rm(2,:));  % Length of vector b (row 2)
    c = norm(Rm(3,:));  % Length of vector c (row 3)
    
    % Store these values in the output structure
    R_struct.a = a;
    R_struct.b = b;
    R_struct.c = c;

    % Precompute dot products to avoid redundant calculations
    dot_b_c = dot(Rm(2,:), Rm(3,:)); % Dot product of b and c
    dot_c_a = dot(Rm(3,:), Rm(1,:)); % Dot product of c and a
    dot_a_b = dot(Rm(1,:), Rm(2,:)); % Dot product of a and b
    
    % Compute angles alpha, beta, gamma in degrees
    R_struct.alpha = acos(dot_b_c / (b * c)) * (180 / pi); % Angle between b and c
    R_struct.beta  = acos(dot_c_a / (c * a)) * (180 / pi); % Angle between c and a
    R_struct.gamma = acos(dot_a_b / (a * b)) * (180 / pi); % Angle between a and b
end
