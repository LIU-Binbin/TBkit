% Define symbolic variables for the components of the wave vector
syms k_x k_y k_z real  % Assuming k_x, k_y, k_z are real

% Momentum combinations in the xy-plane
K_plus  = k_x + 1i* k_y;    % K+ (complex combination of k_x and k_y)
K_minus = k_x - 1i* k_y;    % K- (complex conjugate of K+)
K_plus_2 = (k_x + 1i* k_y)^2;  % Square of K+
K_minus_2 = (k_x - 1i* k_y)^2; % Square of K-

% Terms involving k_x and k_y
K_x2y2 = k_x^2 + k_y^2;  % k_x^2 + k_y^2 (magnitude squared in the xy-plane)

% Higher-order terms
K_pm_3 = 2*k_x*(k_x^2 - 3*k_y^2);  % A higher-order momentum term involving k_x and k_y

% Individual momentum components
K_x = k_x;  % Component k_x
K_y = k_y;  % Component k_y
K_z = k_z;  % Component k_z (not used in this part)

% Squared components
K_x2 = k_x^2;  % k_x^2
K_y2 = k_y^2;  % k_y^2
K_z2 = k_z^2;  % k_z^2 (not used in this part)

% New momentum combinations (likely for a hexagonal or 2D model)
K_1 = k_x;  % k_x term
K_2 = (k_x + sqrt(3)*k_y)/2;  % A rotated momentum component in the xy-plane
K_12 = (k_x - sqrt(3)*k_y)/2; % Another rotated momentum component in the xy-plane
