function [orbital_out, Rm_s_fin] = AddVacuumLayer(orbital_init, POSCAR_file, fin_dir_list, options)
%ADDVACUUMLAYER Add vacuum layer to a structure and adjust orbital coordinates.
%   [ORBITAL_OUT, RM_S_FIN] = ADDVACUUMLAYER(ORBITAL_INIT, POSCAR_FILE, FIN_DIR_LIST, OPTIONS) 
%   modifies lattice vectors to include a vacuum layer and adjusts orbital coordinates 
%   accordingly. This is typically used for surface/interface modeling to avoid periodic 
%   boundary effects.
%
%   Input Arguments:
%   - ORBITAL_INIT: Initial orbital fractional coordinates (Nx3 matrix).
%   - POSCAR_FILE: Path to POSCAR file or direct lattice matrix (3x3). Default: 'POSCAR'.
%   - FIN_DIR_LIST: 1x3 vector specifying vacuum direction weights [a,b,c]. Default: [1 1 0].
%   - OPTIONS: Optional name-value pairs:
%       * fast: Boolean flag for fast mode (skip supercell generation). Default: true.
%       * vacuum_length: Vacuum layer thickness in Angstroms. Default: 10.
%
%   Output:
%   - ORBITAL_OUT: Adjusted orbital coordinates in new fractional basis (Nx3 matrix).
%   - RM_S_FIN: Modified lattice vectors with vacuum layer (3x3 matrix).
%
%   Example Usage:
%   % Basic usage with default parameters
%   [new_orbitals, new_lattice] = AddVacuumLayer(init_orbitals);
%
%   % Custom vacuum direction and length
%   opts = struct('fast', false, 'vacuum_length', 15);
%   [new_orbitals, new_lattice] = AddVacuumLayer(init_orbitals, 'myPOSCAR', [1 0 0], opts);
%
%   Notes:
%   - In fast mode (default), lattice is directly extended without supercell reconstruction.
%   - When POSCAR_FILE is a numeric matrix, it is interpreted as direct lattice vectors.
%   - Coordinate transformation preserves orbital positions in real space.
%
%   See also: POSCAR_READ, SUPERCELL.

% Argument validation block
arguments
    orbital_init
    POSCAR_file = 'POSCAR'
    fin_dir_list = [1 1 0]
    options.fast logical = true
    options.vacuum_length double = 10
end

% Initialize output
orbital_out = orbital_init;
fin_orb = orbital_out;
Ns = eye(3); % Standard basis matrix

% Lattice processing
if ~options.fast
    % Full supercell generation path
    [Rm_tmp, sites_tmp, Atom_name_tmp, Atom_num_tmp] = TBkit.POSCAR_read(POSCAR_file);
    H_hr.supercell(Ns, 'POSCAR_super_fin', Rm_tmp, sites_tmp, Atom_name_tmp, Atom_num_tmp, fin_dir_list);
else
    % Fast mode: direct lattice handling
    switch class(POSCAR_file)
        case {'string', 'char'}
            [Rm_tmp, ~, ~, ~] = TBkit.POSCAR_read(POSCAR_file);
        case 'double'
            Rm_tmp = POSCAR_file;
    end
end

% Modify lattice vectors
Rm_tmp = Ns * Rm_tmp; % Expand lattice
Rmlength = [norm(Rm_tmp(1,:)), norm(Rm_tmp(2,:)), norm(Rm_tmp(3,:))];
Rm_s_fin_add = options.vacuum_length * ...
    [fin_dir_list(1)*Rm_tmp(1,:)/Rmlength(1); ...
     fin_dir_list(2)*Rm_tmp(2,:)/Rmlength(2); ...
     fin_dir_list(3)*Rm_tmp(3,:)/Rmlength(3)];
Rm_s_fin = Rm_tmp + Rm_s_fin_add;

% Coordinate transformation
Rc_add = [0.5, 0.5, 0.5]; % Centroid shift
Rr_add = Rc_add * Rm_s_fin_add;

for i = 1:size(fin_orb, 1)
    Rr_orb = fin_orb(i,:) * Rm_tmp;    % Real-space position
    Rr_new = Rr_orb + Rr_add;          % Apply vacuum shift
    fin_orb(i,:) = Rr_new / Rm_s_fin;  % Convert back to fractional
end

orbital_out = fin_orb;
end