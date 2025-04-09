%% MATLAB Tutorial: Band Structure Plotting with pbandplot()

%% Basic Band Plot Configuration
% Plot bands within [-6,6] energy window
% nosoc
pbandplot('Ecut',[-6,6]);
% 'Ecut' specifies energy range for visualization
% Default plotting style will be used

%% PROCAR Data Processing
% Read PROCAR file with spin-orbit coupling (SOC) enabled
[PROCAR_collection,EIGENCAR] = PROCAR_read('SOC_flag',1);
% PROCAR_collection: Cell array of orbital-projected weights
% EIGENCAR: Eigenvalues matrix
% SOC_flag=1 enables spin-orbit coupling processing

% Set Fermi energy (system-specific value)
Fermi = -2.97452948; % Adjust according to your calculation

%% Orbital Weight Configuration
% Define orbital weights for visualization
% Column indices correspond to:
% s  py pz px dxy dyz dz2 dxz x2-y2 fy3x2 fxyz fyz2 fz3 fxz2 fzx2 fx3 tot
WEIGHTCAR_C1(1).WEIGHTCAR = PROCAR_collection{1,1}; % s-orbital (Column 1)
WEIGHTCAR_C1(2).WEIGHTCAR = PROCAR_collection{1,3}; % pz-orbital (Column 3)

% Set display labels
WEIGHTCAR_C1(1).displayname = 'buckled Graphene:C_1-s';
WEIGHTCAR_C1(2).displayname = 'buckled Graphene:C_1-p_z';

%% Advanced Band Plotting
% Generate customized band structure plot
pbandplot(WEIGHTCAR_C1, EIGENCAR,...
    'silent', true,...    % Suppress command window output
    'Ecut',[-15,15]);     % Extended energy range visualization

% Key Features:
% 1. Orbital-projected weights shown as color maps
% 2. Multiple orbitals plotted simultaneously
% 3. Fermi level alignment based on input value

% WEIGHTCAR datatype
% ************************************************************************************************************************************************
%     ----- band     kpoint1     kpoint2     kpoint3     ............
%           band1                                .
%           band2                                .
%             .       .     .      .      .   weight      .     .     .     .
%             .                                  .
%             .                                  .
% ************************************************************************************************************************************************


