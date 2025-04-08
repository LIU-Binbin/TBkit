%% CuI Projected Band Structure Tutorial
% This tutorial demonstrates how to plot projected band structure for Copper Iodide (CuI)

%% Step 1: Data Preparation
% Ensure these files/folders exist in current path:
% - PBAND*.dat      : VASPKIT-generated pband weight data
% - EIGENVAL        : Energy eigenvalues file
% - KPOINTS
% - POSCAR
% - DOSCAR

%% Step 2: Read Weight Data
[WEIGHTCAR_struct_cell, Name_list, ~] = WEIGHTCAR_read_dir('PBAND');

%% Step 3: Atomic Orbital Processing
% Separate atomic orbital data
WEIGHTCAR_struct_Cu = WEIGHTCAR_struct_cell{1};  % Copper orbitals
WEIGHTCAR_struct_I1 = WEIGHTCAR_struct_cell{2};  % Iodine-1 orbitals
WEIGHTCAR_struct_I2 = WEIGHTCAR_struct_cell{3};  % Iodine-2 orbitals

%% Step 4: Orbital Configuration
% First configuration scheme
WEIGHTCAR_struct_plot(1) = WEIGHTCAR_struct_Cu(1);      % Cu-d orbitals
WEIGHTCAR_struct_plot(2) = WEIGHTCAR_struct_I1(2);      % I1-py
WEIGHTCAR_struct_plot(3) = WEIGHTCAR_struct_I1(4);      % I1-px

% Second configuration (px+py combined)
WEIGHTCAR_struct_plot_2(1) = WEIGHTCAR_struct_Cu(1);    % Cu-d
WEIGHTCAR_struct_plot_2(2) = WEIGHTCAR_struct_I1(2);    % I1-p
WEIGHTCAR_struct_plot_2(2).displayname = 'I1-px,py';    % Display name
WEIGHTCAR_struct_plot_2(2).WEIGHTCAR = ...             % Combine px & py
    WEIGHTCAR_struct_I1(2).WEIGHTCAR + ...             
    WEIGHTCAR_struct_I1(4).WEIGHTCAR;
WEIGHTCAR_struct_plot_2(2).WEIGHTCAR = ...             % Normalization
    WEIGHTCAR_struct_plot_2(2).WEIGHTCAR / 2;

%% Step 5: Band Structure Visualization
EIGENCAR = EIGENVAL_read;       % Read eigenvalues
Ecut = [-3, 3];                 % Energy range (eV)
titlestring = 'fatband-CuI-nosoc'; 

% Generate band plot
[fig, ax] = pbandplot(WEIGHTCAR_struct_plot_2, EIGENCAR,...
    'Ecut', Ecut,...
    'title', titlestring,...
    'silent', true);

%% Parameter Reference
% -------------------------------
% | Parameter    | Description          |
% |--------------|----------------------|
% | WEIGHTCAR    | Orbital weight matrix (kpoints×bands) |
% | EIGENCAR     | Energy eigenvalues matrix |
% | Ecut         | Energy display range [min, max] |
% | silent       | Disable command line output |
% -------------------------------
