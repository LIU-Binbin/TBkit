function dirstring = nanodisk_vasp_prework(level_cut, search_range)
    % NANODISK_VASP_PREWORK Preprocesses nanodisk system for VASP calculations.
    %
    %   dirstring = NANODISK_VASP_PREWORK(level_cut, search_range) performs
    %   preprocessing tasks for a nanodisk system, including reading the POSCAR
    %   file, generating a tight-binding Hamiltonian, and saving the results.
    %
    %   Inputs:
    %       level_cut    - Cutoff level for the tight-binding Hamiltonian.
    %       search_range - Range for nearest neighbor search.
    %
    %   Outputs:
    %       dirstring    - Name of the saved MAT-file containing the results.

    % Validate and set default values for input arguments
    if nargin < 1
        level_cut = 2; % Default cutoff level
    end
    if nargin < 2
        search_range = [0, 0, 0]; % Default search range
    end

    % Read the POSCAR file to obtain lattice vectors and atomic positions
    POSCAR_read;

    % Perform nearest neighbor search to obtain atomic neighbors
    [Atom_store, nn_store, Rnn] = nn(Rm, sites, search_range);

    % Generate the tight-binding Hamiltonian for the nanodisk system
    H = H_TB_gen(level_cut, nn_store, sites, 'nano');

    % Save the results to a MAT-file
    current_dir = pwd;
    dir_parts = strsplit(current_dir, '/');
    dirstring = strcat(dir_parts{end}, '.mat');
    save(dirstring);
end
