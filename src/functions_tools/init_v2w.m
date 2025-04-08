function dirstring = init_v2w(SOCFLAG, Erange, options)
% init_v2w - Initialize and prepare files for Wannier90 calculations.
%
% Usage:
%   dirstring = init_v2w(SOCFLAG, Erange, options)
%   dirstring = init_v2w(SOCFLAG, Erange)
%   dirstring = init_v2w(SOCFLAG)
%   dirstring = init_v2w()
%
% Inputs:
%   SOCFLAG      - Flag for Spin-Orbit Coupling (SOC), default is 0 (no SOC).
%   Erange       - Energy range for the bands, default is [-3, 3].
%   options      - A struct with additional options:
%     - check        : If true, prompts user to check and find the bands.
%     - Projector_list: A list of projectors for the Wannier90 calculation.
%
% Outputs:
%   dirstring    - A string representing the directory path.

    %% Set default arguments if not provided
    arguments
        SOCFLAG = 0;  % Default to no Spin-Orbit Coupling (SOC)
        Erange = [-3, 3];  % Default energy range for bands
        options.check = false;  % Default to no check
        options.Projector_list = [];  % Default empty projector list
    end

    %% Import custom functions (if necessary)
    import park.*;

    %% Read POSCAR file for lattice and atom information
    SOCflag = SOCFLAG;
    [Rm, sites, Atom_name, Atom_num] = POSCAR_read();

    %% Try reading EIGENCAR for band structure information
    try
        check = true;
        EIGENCAR = EIGENVAL_read();  % Read band eigenvalues from EIGENVAL file
    catch
        EIGENCAR = [];  % If EIGENCAR is not found, set it to an empty array
        check = false;  % Set check flag to false if EIGENCAR is not available
    end

    %% Optionally check and find bands
    if options.check && nargin < 2
        while ~strcmp(input('y/n for continue or break to find bands', 's'), 'n')
            findbands(EIGENCAR);  % Function to find bands
        end
    elseif options.check
        findbands(EIGENCAR, Erange);  % If check is true, find bands within the given energy range
    end

    %% Generate kpath for band structure visualization
    kpath_card_gen();

    %% Write projection information
    if isempty(options.Projector_list)
        num_wan = write_pj();  % If no projector list provided, use default projector
    else
        num_wan = write_pj(options.Projector_list);  % Use provided projector list
    end

    %% Generate the Wannier90 input files with SOC and band selection
    wt_in_gen(Rm, sites, Atom_num, Atom_name, SOCflag);  % Generate the Wannier90 input for SOC

    % If the EIGENCAR file was available and the check flag is set, select bands
    if check && options.check
        [num_bands, Nmin, Nmax] = selectbands(EIGENCAR, Rm);  % Select bands for Wannier90
        Nbands = size(EIGENCAR, 1);  % Get the total number of bands
        wannier90_win_gen(SOCflag, num_bands, num_wan, Nmin, Nmax, Nbands);  % Generate Wannier90 input with selected bands
    else
        wannier90_win_gen(SOCflag, 10, num_wan, 1, 10, 10);  % Generate default Wannier90 input if no bands are selected
    end

    %% Completion message
    disp('All done');

end
