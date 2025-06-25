function [pgon, vertices] = generateKochSnowflake(n, ~, center, radius)
%GENERATEKOCHSNOWFLAKE  Create a Sierpinski triangle polyshape and its vertices
%
%   [PGON, VERTICES] = GENERATEKOCHSNOWFLAKE(N, ~, CENTER, RADIUS)
%   generates a Sierpinski triangle by recursively subdividing an equilateral
%   triangle into smaller triangles and removing the central one.
%
%   INPUTS
%     N       – (scalar) Number of recursion levels (non-negative integer).
%               Default: 3.
%     THETA   – (ignored) Retained for interface compatibility.
%     CENTER  – (1×2)   [cx, cy] coordinates of the center.
%               Default: [0, 0].
%     RADIUS  – (scalar) Circumradius of the initial equilateral triangle.
%               Default: 10.
%
%   OUTPUTS
%     PGON     – polyshape object containing the filled Sierpinski triangle.
%     VERTICES – (M×2) array of vertex coordinates in drawing order (with NaN separators).
%
%   Example
%     [pg, vtx] = generateKochSnowflake(4);   % 4 recursion levels
%     plot(pg); axis equal off
%
%   See also POLYSHAPE, PATCH

    % ---------- Default arguments ---------------------------------------
    if nargin < 1 || isempty(n),      n      = 3;      end
    % theta is ignored (retained for interface compatibility)
    if nargin < 3 || isempty(center), center = [0, 0]; end
    if nargin < 4 || isempty(radius), radius = 10;     end

    % ---------- Build the base equilateral triangle ---------------------
    baseAngles = (0:2).' * 2*pi/3 + pi/2;           % 90°, 210°, 330°
    baseTri    = radius * [cos(baseAngles), ...
                           sin(baseAngles)] + center;
    
    % ---------- Recursive construction of Sierpinski triangle -----------
    % Initialize with the base triangle
    triangles = {baseTri};  % Cell array to store all triangles
    
    % Perform recursive subdivision
    for k = 1:n
        newTriangles = cell(1, 3*numel(triangles)); % Each triangle produces 3
        triCount = 0;
        
        for idx = 1:numel(triangles)
            % Get vertices of current triangle
            tri = triangles{idx};
            A = tri(1,:);
            B = tri(2,:);
            C = tri(3,:);
            
            % Calculate midpoints of sides
            M_AB = (A + B)/2;
            M_BC = (B + C)/2;
            M_CA = (C + A)/2;
            
            % Create three new sub-triangles
            newTriangles{triCount+1} = [A; M_AB; M_CA];   % Top sub-triangle
            newTriangles{triCount+2} = [M_AB; B; M_BC];   % Right sub-triangle
            newTriangles{triCount+3} = [M_CA; M_BC; C];   % Left sub-triangle
            
            triCount = triCount + 3;
        end
        triangles = newTriangles(1:triCount);
    end
    
    % ---------- Convert to polyshape format ----------------------------
    % Preallocate vertex array with NaN separators
    numTriangles = numel(triangles);
    vertices = NaN(4*numTriangles, 2);  % 3 vertices + NaN per triangle
    
    % Populate vertex array
    for idx = 1:numTriangles
        startIdx = (idx-1)*4 + 1;
        tri = triangles{idx};
        
        % Add triangle vertices (closed polygon)
        vertices(startIdx:startIdx+2, :) = tri;
        vertices(startIdx+3, :) = [NaN, NaN];  % Separator
    end
    
    % Remove trailing NaN if exists
    if all(isnan(vertices(end,:)))
        vertices = vertices(1:end-1,:);
    end
    
    % Create polyshape object
    pgon = polyshape(vertices(:,1), vertices(:,2), 'Simplify', false);
end