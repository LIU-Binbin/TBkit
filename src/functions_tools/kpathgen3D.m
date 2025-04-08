function [klist_cart, klist_l, klist_frac, kpoints_l, kpoints_name] = kpathgen3D(Rm, kpoints, nodes, kpoints_name)   
    % K-Point Path Generation for 3D and 2D Lattice Structures
    % This function generates k-points in Cartesian (real) space, fractional space,
    % linear space, and calculates the path lengths based on input parameters.
    %
    % Inputs:
    %    Rm           - Lattice vectors (real space)
    %    kpoints      - List of k-points (each row is a k-point)
    %    nodes        - Number of nodes for interpolation between k-points
    %    kpoints_name - Names of the k-points (optional)
    %
    % Outputs:
    %    klist_cart   - K-points in Cartesian (real) space
    %    klist_l      - K-points in linear space
    %    klist_frac   - K-points in fractional space (interpolated k-points)
    %    kpoints_l    - List of k-point path lengths
    %    kpoints_name - Names of the k-points (optional)
    
    %% Handle Input Arguments
    % Read POSCAR or KPOINTS file if no inputs are provided
    if nargin < 1
        POSCAR_read;  % Read POSCAR if no arguments are given
    end
    
    if nargin < 2
        [kpoints, nodes, kpoints_name] = KPOINTS_read();  % Read from KPOINTS file if kpoints are not provided
    elseif nargin < 3
        [kpoints, nodes, kpoints_name] = KPOINTS_read(kpoints);  % Read nodes if kpoints are given
    end
    % bug
    dim = max(size(Rm));
    clear length
    %
    %% Initialize Parameters
    % Calculate the reciprocal lattice matrix (Gk)
    Gk = (eye(dim) * 2 * pi / Rm)';  
    n = size(kpoints, 1) / 2;  % Number of k-point segments
    
    %% Handle Nodes (if only one node count is provided)
    if max(size(nodes)) == 1
        nodes = ones(n, 1) * nodes;  % If only one node value is provided, expand it for each segment
    end
    
    %% Generate K-points (Fractional, Linear, and Cartesian)
    if dim == 3  % 3D case (real space)
        klist_frac = [];  % Initialize fractional k-points
        
        % Interpolate k-points in fractional space (3D)
        for i = 1:n
            klisttempX = linspace(kpoints(2*i-1, 1), kpoints(2*i, 1), nodes(i));
            klisttempY = linspace(kpoints(2*i-1, 2), kpoints(2*i, 2), nodes(i));
            klisttempZ = linspace(kpoints(2*i-1, 3), kpoints(2*i, 3), nodes(i));
            klist_frac = [klist_frac; klisttempX' klisttempY' klisttempZ'];  % Append interpolated k-points
        end
        
        % Calculate the path lengths in fractional space
        kpoints_l = [0];  % Initialize with 0 for the first point
        for i = 1:n
            length = norm((kpoints(2*i, :) - kpoints(2*i-1, :)) * Gk);  % Calculate distance in reciprocal space
            kpoints_l = [kpoints_l; kpoints_l(i) + length];  % Update cumulative path length
        end
        
        % Interpolate k-points in linear space
        klist_l = [];
        for i = 1:n
            klisttemp = linspace(kpoints_l(i), kpoints_l(i+1), nodes(i));  % Linearly interpolate
            klist_l = [klist_l klisttemp];  % Append the results
        end
        
        % Calculate k-points in Cartesian (real) space for calculation purposes
        klist_cart = klist_frac * Gk;  % Convert fractional k-points to Cartesian space
        
    elseif dim == 2  % 2D case (real space)
        klist_frac = [];  % Initialize fractional k-points
        
        % Interpolate k-points in fractional space (2D)
        for i = 1:n
            klisttempX = linspace(kpoints(2*i-1, 1), kpoints(2*i, 1), nodes(i));
            klisttempY = linspace(kpoints(2*i-1, 2), kpoints(2*i, 2), nodes(i));
            klist_frac = [klist_frac; klisttempX' klisttempY'];  % Append interpolated 2D k-points
        end
        
        % Calculate the path lengths in fractional space (2D)
        kpoints_l = [0];  % Initialize with 0 for the first point
        for i = 1:n
            length = norm((kpoints(2*i, 1:2) - kpoints(2*i-1, 1:2)) * Gk);  % 2D length calculation
            kpoints_l = [kpoints_l; kpoints_l(i) + length];  % Update cumulative path length
        end
        
        % Interpolate k-points in linear space (2D)
        klist_l = [];
        for i = 1:n
            klisttemp = linspace(kpoints_l(i), kpoints_l(i+1), nodes(i));  % Linearly interpolate
            klist_l = [klist_l klisttemp];  % Append the results
        end
        
        % Calculate k-points in Cartesian (real) space for calculation purposes (2D)
        klist_cart = klist_frac * Gk;  % Convert fractional k-points to Cartesian space
    end
end
