function [klist_cart,klist_frac,klist_l,kpoints_l,kpoints_frac] = kpathgen(kpoints,nodes,Gk,Gk_,options)
%KPATHGEN Generate k-point paths for Brillouin zone sampling
%   Generates both fractional and Cartesian coordinates for k-point paths
%   between specified high-symmetry points, with optional reciprocal lattice conversions
%
%   Inputs:
%   kpoints    : Nx3 matrix of high-symmetry points in fractional coordinates (pairs define path segments)
%   nodes      : Number of points per segment (scalar) or array specifying points per segment
%   Gk         : Reciprocal lattice vectors matrix for Cartesian conversion (3x3)
%   Gk_        : Inverse reciprocal lattice matrix for reverse conversion
%   options.Dim: Dimensionality of the system (default=3)
%
%   Outputs:
%   klist_cart   : Cartesian coordinates of k-points
%   klist_frac   : Fractional coordinates of k-points
%   klist_l      : Linearized path length coordinates
%   kpoints_l    : Cumulative path lengths at segment boundaries
%   kpoints_frac : Original input kpoints (fractional coordinates)

    arguments
        kpoints;
        nodes = 60;
        Gk = [];
        Gk_ = [];
        options.Dim = 3;
    end
    
    % Determine operation mode based on lattice matrices
    if nargin > 3 && ~isequal(Gk,Gk_)
        mode = 'conversion_mode';
    else
        mode = 'normal_mode';
    end
    
    % Initialize system dimensionality
    Dim = options.Dim;
    
    % Validate input dimensions
    if mod(size(kpoints,1),2) ~= 0
        error('kpoints must contain even number of rows (segment pairs)');
    end
    
    % Calculate number of path segments
    num_segments = size(kpoints,1)/2;
    
    % Initialize nodes array
    if isscalar(nodes)
        nodes = repmat(nodes, num_segments, 1);
    end
    
    % Preallocate fractional coordinate array
    klist_frac = zeros(sum(nodes), Dim);
    
    % Generate k-point paths for each segment
    for seg_idx = 1:num_segments
        % Extract start and end points for current segment
        start_point = kpoints(2*seg_idx-1, :);
        end_point = kpoints(2*seg_idx, :);
        
        % Create linearly spaced points for each dimension
        seg_points = arrayfun(@(d) linspace(start_point(d), end_point(d), nodes(seg_idx))',...
            1:Dim, 'UniformOutput', false);
        
        % Store generated points
        start_idx = sum(nodes(1:seg_idx-1)) + 1;
        end_idx = sum(nodes(1:seg_idx));
        klist_frac(start_idx:end_idx, :) = [seg_points{:}];
    end
    
    % Convert to Cartesian coordinates using reciprocal lattice vectors
    klist_cart = klist_frac * Gk;
    
    % Calculate cumulative path lengths for visualization
    kpoints_l = zeros(num_segments+1, 1);
    for seg_idx = 1:num_segments
        segment_vector = (kpoints(2*seg_idx,:) - kpoints(2*seg_idx-1,:)) * Gk;
        kpoints_l(seg_idx+1) = kpoints_l(seg_idx) + norm(segment_vector);
    end
    
    % Generate linearized path coordinates
    klist_l = zeros(1, sum(nodes));
    for seg_idx = 1:num_segments
        seg_lin = linspace(kpoints_l(seg_idx), kpoints_l(seg_idx+1), nodes(seg_idx));
        start_idx = sum(nodes(1:seg_idx-1)) + 1;
        end_idx = sum(nodes(1:seg_idx));
        klist_l(start_idx:end_idx) = seg_lin;
    end
    
    % Handle coordinate system conversion if required
    if strcmp(mode, 'conversion_mode')
        klist_frac = klist_cart / Gk_;
    end
    
    % Preserve original fractional points
    kpoints_frac = kpoints;
end