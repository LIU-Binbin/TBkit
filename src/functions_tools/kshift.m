function [klist_new, DATA_L_new] = kshift(klist, KCUBE_BULK, DATA_L, opt)
% KSHIFT shifts or expands the k-point mesh and interpolates data at new k-points.
%   Input:
%     klist        : Array of k-points (N x 3).
%     KCUBE_BULK   : Bulk lattice vectors (N x 3).
%     DATA_L       : Data associated with k-points (N x M).
%     opt          : Structure of optional parameters, including:
%       cart      : Boolean, whether to convert k-points to fractional coordinates.
%       Rm        : Lattice matrix (reciprocal).
%       tol       : Tolerance for uniqueness of k-points (default: 1e-12).
%       MeshOutput: Boolean, whether to output mesh data (default: false).
%       method    : Interpolation method (options: 'linear', 'nearest', 'natural').
%   Output:
%     klist_new   : Expanded list of k-points (new k-points).
%     DATA_L_new  : Interpolated data at the new k-points.

% Default arguments
arguments
    klist        = [];  % K-point list (default empty)
    KCUBE_BULK   = [-0.5 -0.5 -0.5; 1 0 0; 0 1 0; 0 0 1];  % Default lattice vectors
    DATA_L       = [];  % Associated data (default empty)
    opt.cart      = false;  % Whether to use Cartesian coordinates (default false)
    opt.Rm        = POSCAR_read();  % Reciprocal lattice matrix (from POSCAR by default)
    opt.tol       = 1e-12;  % Tolerance for uniqueness of k-points
    opt.MeshOutput = false;  % Mesh output flag
    opt.method    {mustBeMember(opt.method, {'linear', 'nearest', 'natural'})} = 'linear';  % Interpolation method
end

%% Reshape DATA_L if needed
sizemesh = size(DATA_L);
if numel(sizemesh) == 3
    DATA_L = reshape(DATA_L, [sizemesh(1) * sizemesh(2), sizemesh(3)]);
elseif (size(klist, 1) * size(klist, 2) == sizemesh(1) * sizemesh(2)) || (size(klist, 1) == sizemesh(1) * sizemesh(2))
    DATA_L = reshape(DATA_L, [sizemesh(1) * sizemesh(2), 1]);
end

% Handle 3D klist (reshape if needed)
sizemeshk = size(klist);
if numel(sizemeshk) == 3
    klist = reshape(klist, [sizemesh(1) * sizemesh(2), 3]);
end

%% Convert klist to fractional coordinates if needed
if opt.cart
    Gk = (2 * pi * eye(3) / opt.Rm).';  % Reciprocal lattice matrix for Cartesian coordinates
    klist = klist / Gk;  % Convert to fractional coordinates
else
    Gk = (2 * pi * eye(3) / opt.Rm).';  % Reciprocal lattice matrix (for fractional coordinates)
end

%% Remove duplicate k-points
[klist, uniqueseq, ~] = uniquetol(klist, opt.tol, 'ByRows', true);
if ~isempty(DATA_L)
    DATA_L = DATA_L(uniqueseq, :);
end

%% Adjust KCUBE_BULK based on lattice dimension (2D or 3D)
switch size(KCUBE_BULK, 1)
    case 2  % 2D lattice (treat as 2D)
        KCUBE_BULK(3, :) = KCUBE_BULK(2, :) + KCUBE_BULK(1, :);  % Generate the third vector
        CENTER = mean(KCUBE_BULK, 1);  % Center of the bulk
        CENTER_base = round(CENTER);  % Shift to the origin
        KCUBE_BULK = KCUBE_BULK(:, 1:2);  % Only keep 2D components
        mode = '2D';
        
    case 3  % 3D lattice
        KCUBE_BULK(4, :) = KCUBE_BULK(2, :) + KCUBE_BULK(3, :);
        KCUBE_BULK(5, :) = KCUBE_BULK(2, :) + KCUBE_BULK(4, :);
        KCUBE_BULK(6, :) = KCUBE_BULK(3, :) + KCUBE_BULK(4, :);
        KCUBE_BULK(7, :) = KCUBE_BULK(2, :) + KCUBE_BULK(3, :) + KCUBE_BULK(4, :);
        CENTER = mean(KCUBE_BULK([1, 7], :), 1);  % Center of the bulk
        CENTER_base = round(CENTER);  % Shift to the origin
        mode = '3D';
end

%% Create a polyhedron (convex hull of KCUBE_BULK)
DT = delaunayTriangulation(KCUBE_BULK);  % Delaunay triangulation for polyhedron
[~, v] = convexHull(DT);  % Compute the convex hull

%% Generate new k-points based on the convex hull
Vceil = ceil(abs(v));  % Max integer range for tiling
VceilL = (-Vceil:Vceil).';  % Create the range of multiples
nVceilL = numel(VceilL);  % Number of points in the expanded mesh

switch mode
    case '2D'  % Expand in 2D
        SuperL = zeros(nVceilL^2, 3);
        count = 0;
        for i = -Vceil:Vceil
            for j = -Vceil:Vceil
                count = count + 1;
                SuperL(count, :) = [i, j, 0] - CENTER_base;  % Offset the lattice center
            end
        end
        
        % Expand the k-list by tiling
        nklist = size(klist, 1);
        SuperL = kron(SuperL, ones([nklist, 1]));  % Replicate for each k-point
        klist_new = repmat(klist, [count, 1]) + SuperL;  % Add shifted k-points
        
        % Find the points inside the convex hull
        ID = pointLocation(DT, klist_new(:, 1:2));
        
    case '3D'  % Expand in 3D
        SuperL = zeros(nVceilL^3, 3);
        count = 0;
        for i = -Vceil:Vceil
            for j = -Vceil:Vceil
                for k = -Vceil:Vceil
                    count = count + 1;
                    SuperL(count, :) = [i, j, k] - CENTER_base;
                end
            end
        end
        
        % Expand the k-list by tiling
        nklist = size(klist, 1);
        SuperL = kron(SuperL, ones([nklist, 1]));  % Replicate for each k-point
        klist_new = repmat(klist, [count, 1]) + SuperL;  % Add shifted k-points
        
        % Find the points inside the convex hull
        ID = pointLocation(DT, klist_new);
end

%% Filter valid k-points and map them back to [0, 1)
klist_new = klist_new(~isnan(ID), :);  % Remove points outside the convex hull
klist_new_map = mod(klist_new, 1);  % Map to unit cell

%% Interpolate data at new k-points (if provided)
if ~isempty(DATA_L)
    [~, reseqL] = ismembertol(klist_new_map, klist, opt.tol, 'ByRows', true);
    rmlist = reseqL == 0;  % Identify points to be removed
    
    if ~sum(rmlist)  % If no points to remove
        DATA_L_new = DATA_L(reseqL, :);
    else
        if ~opt.MeshOutput
            switch mode
                case '2D'  % Interpolate in 2D
                    DATA_L_new = zeros(size(klist_new_map, 1), size(DATA_L, 2));
                    klist1 = klist(:, 1); klist2 = klist(:, 2);
                    klistnew1 = klist_new_map(:, 1); klistnew2 = klist_new_map(:, 2);
                    for i = 1:size(DATA_L, 2)
                        F = scatteredInterpolant(klist1, klist2, DATA_L(:, i), opt.method);
                        DATA_L_new(:, i) = F(klistnew1, klistnew2);
                    end
                case '3D'  % Interpolate in 3D
                    DATA_L_new = zeros(size(klist_new_map, 1), size(DATA_L, 2));
                    klist1 = klist(:, 1); klist2 = klist(:, 2); klist3 = klist(:, 3);
                    klistnew1 = klist_new_map(:, 1); klistnew2 = klist_new_map(:, 2); klistnew3 = klist_new_map(:, 3);
                    for i = 1:size(DATA_L, 2)
                        F = scatteredInterpolant(klist1, klist2, klist3, DATA_L(:, i));
                        DATA_L_new(:, i) = F(klistnew1, klistnew2, klistnew3);
                    end
            end
        end
    end
end

%% Convert klist_new back to Cartesian coordinates if needed
if opt.cart
    klist_new = klist_new * Gk;  % Convert back to Cartesian
end

end
