function [pgon, vertices] = generateKochSnowflake(n, theta, center, radius)
%GENERATEKOCHSNOWFLAKE  Create a Koch–snowflake polyshape and its vertices
%
%   [PGON, VERTICES] = GENERATEKOCHSNOWFLAKE(N, THETA, CENTER, RADIUS)
%   generates a Koch snowflake by recursively replacing every line segment
%   with four segments that form an outward-pointing isosceles triangle.
%
%   INPUTS
%     N       – (scalar) Number of recursion levels (non-negative integer).
%               Default: 3.
%     THETA   – (scalar) Interior rotation angle in radians used to create
%               the "bump" at each recursion step.  THETA = pi/3 (60°)
%               reproduces the classical Koch snowflake, which exhibits
%               six-fold (hexagonal) rotational symmetry.  Default: pi/3.
%     CENTER  – (1×2)   [cx, cy] coordinates of the snowflake centre.
%               Default: [0, 0].
%     RADIUS  – (scalar) Circumradius of the initial equilateral triangle.
%               Default: 10.
%
%   OUTPUTS
%     PGON     – polyshape object containing the filled snowflake.
%     VERTICES – (M×2) array of vertex coordinates in drawing order.
%
%   Example
%     [pg, vtx] = generateKochSnowflake(4);   % 4 recursion levels
%     plot(pg); axis equal off
%
%   See also POLYSHAPE, PATCH
% -------------------------------------------------------------------------

    % ---------- Default arguments ---------------------------------------
    if nargin < 1 || isempty(n),      n      = 3;      end
    if nargin < 2 || isempty(theta),  theta  = -pi/3;   end   % 60° for hexagon
    if nargin < 3 || isempty(center), center = [0, 0]; end
    if nargin < 4 || isempty(radius), radius = 10;     end

    % ---------- Build the base equilateral triangle ---------------------
    baseAngles   = (0:2).' * 2*pi/3 + pi/2;           % 90°, 210°, 330°
    baseTri      = radius * [cos(baseAngles), ...
                              sin(baseAngles)] + center;
    x = baseTri(:,1);                                 % x-coordinates
    y = baseTri(:,2);                                 % y-coordinates
    x(end+1) = x(1);                                  % close polygon
    y(end+1) = y(1);

    % ---------- Recursive construction ----------------------------------
    for k = 1:n
        [x, y] = kochIteration(x, y, theta);
    end

    vertices = [x(:), y(:)];
    pgon     = polyshape(vertices);                   % build polyshape
end
% ======================================================================= %
function [xNew, yNew] = kochIteration(x, y, theta)
%KOCHITERATION  Perform one Koch subdivision on a vertex list.

    xNew = [];                                       % preallocate
    yNew = [];

    rotMat = [cos(theta), -sin(theta);               % rotation matrix
              sin(theta),  cos(theta)];

    for k = 1:numel(x)-1
        % Endpoints of current segment
        x0 = x(k);   y0 = y(k);
        x1 = x(k+1); y1 = y(k+1);

        % Vector from P0 to P1
        dx = x1 - x0;
        dy = y1 - y0;

        % Points at 1/3 and 2/3 along the segment
        xA = x0 + dx/3;        yA = y0 + dy/3;
        xB = x0 + 2*dx/3;      yB = y0 + 2*dy/3;

        % Outward-pointing vertex (rotate the middle third)
        vRot = rotMat * ([xB; yB] - [xA; yA]);
        xC   = xA + vRot(1);
        yC   = yA + vRot(2);

        % Append four points that replace the old segment
        xNew = [xNew; x0; xA; xC; xB];
        yNew = [yNew; y0; yA; yC; yB];
    end

    % Close the polyline
    xNew = [xNew; x(end)];
    yNew = [yNew; y(end)];
end
