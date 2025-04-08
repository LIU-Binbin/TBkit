function maprule = orbital_maprule_gen(f_or_not)
    % orbital_maprule_gen generates a mapping of integers to orbital types.
    %
    % Inputs:
    %   f_or_not - Optional flag to determine the range of orbital types.
    %              If 1, maps integers 1-17 to orbital types 's' to 'fx3'.
    %              If 0 or not provided, maps integers 1-10 to orbital types 's' to 'dx2-y2'.
    %
    % Outputs:
    %   maprule  - A containers.Map object mapping integers to orbital types.
    
    % Set default value for f_or_not if not provided
    if nargin < 1
        f_or_not = 0;
    end
    
    % Define orbital types based on f_or_not
    if f_or_not == 1
        % Map integers 1-17 to orbital types 's' to 'fx3'
        orbital_types = {'s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2', ...
                         'fy3x2', 'fxyz', 'fyz2', 'fz3', 'fxz2', 'fzx2', 'fx3', 'tot'};
        keys = 1:17;
    else
        % Map integers 1-10 to orbital types 's' to 'dx2-y2'
        orbital_types = {'s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2', 'tot'};
        keys = 1:10;
    end
    
    % Create the containers.Map object
    maprule = containers.Map(keys, orbital_types);
end
