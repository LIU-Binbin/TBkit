function maprule = orbital_maprule_gen(f_or_not)
    % ORBITAL_MAPRULE_GEN Generates a mapping of orbital types based on the input.
    %
    % Usage:
    %   maprule = orbital_maprule_gen(f_or_not)
    %
    % Inputs:
    %   f_or_not - (Optional) A flag indicating the type of mapping:
    %              0 (default) - maps first 10 orbitals
    %              1           - maps 17 orbitals
    %
    % Outputs:
    %   maprule - A containers.Map object mapping integers to orbital types.

    % Check if the input argument is provided, defaulting to 0 if not
    if nargin < 1
        f_or_not = 0;
    end
    
    % Create the mapping based on the input value
    if f_or_not == 1
        % Create a map for a complete set of orbitals (1 to 17)
        maprule = containers.Map( ...
            {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17}, ...
            {'s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2', ...
             'fy3x2', 'fxyz', 'fyz2', 'fz3', 'fxz2', 'fzx2', 'fx3', 'tot'});
    elseif f_or_not == 0
        % Create a map for a subset of orbitals (1 to 10)
        maprule = containers.Map( ...
            {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, ...
            {'s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2', 'tot'});
    end
end