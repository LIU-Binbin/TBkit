function varargout = BCplot2D(BCCAR, Grid, Rm, options)
    % BCplot2D - A function to visualize 2D band structure or related data.
    %
    % Usage:
    %   BCplot2D(BCCAR, Grid, Rm, options)
    %
    % Inputs:
    %   BCCAR   - 2D or 3D matrix containing the data to be plotted.
    %   Grid    - Grid points (coordinates) associated with BCCAR data.
    %   Rm      - Lattice constant (used for plotting reciprocal space).
    %   options - A structure containing optional parameters for customization:
    %     - ColorCut (default 1): Scale factor for the maximum value in the color map.
    %     - ColorCutMinus (default -1): Scale factor for the minimum value in the color map.
    %     - ax (default empty): Handle to a specific axes object for plotting.
    %     - BZ (default false): Whether to plot Brillouin Zone.
    %     - BZmode (default '3D'): Mode for the Brillouin Zone ('3D' or '2D').
    %     - BZlabel (default true): Whether to label the Brillouin Zone.
    %     - shading (default false): Whether to apply shading to the plot.
    %     - scalefactor (default 100): Scaling factor for data visualization.
    %     - GkOrigin (default [0 0 0]): Origin point for the reciprocal lattice.
    %     - method (default 'linear'): Interpolation method for meshgrid (options: 'linear', 'nearest', 'natural', 'v4').
    %     - plotmode (default 'surf'): Plotting mode ('surf' or 'contour').
    %
    % Outputs:
    %   varargout - Depending on the number of output arguments requested, 
    %              this can return the plot axes or handle objects.

    arguments
        BCCAR double = [];                    % BCCAR data matrix.
        Grid double = [];                     % Grid points.
        Rm double = [];                       % Lattice constant for reciprocal space.
        options.ColorCut double = 1;          % Color scaling factor for positive values.
        options.ColorCutMinus double = -1;    % Color scaling factor for negative values.
        options.ax handle = handle([]);       % Axes handle for plotting.
        options.BZ logical = false;           % Flag to plot Brillouin Zone.
        options.BZmode {mustBeMember(options.BZmode, {'3D', '2D'})} = '3D';  % Brillouin Zone mode ('3D' or '2D').
        options.BZlabel logical = true;       % Flag to add Brillouin Zone labels.
        options.shading logical = false;      % Whether to apply shading to the plot.
        options.scalefactor double = 100;     % Scaling factor for data visualization.
        options.GkOrigin double = [0, 0, 0]; % Origin point for the reciprocal lattice.
        options.method {mustBeMember(options.method, {'linear', 'nearest', 'natural', 'v4'})} = 'linear'; % Interpolation method.
        options.plotmode {mustBeMember(options.plotmode, {'surf', 'contour'})} = 'surf'; % Plotting mode ('surf' or 'contour').
    end

    % Handle missing axes input and create a new figure if necessary
    if isempty(options.ax)
        Fig = Figs(1, 1);  % Create a new figure window
        ax = Fig.axes(1);  % Use the first axis of the new figure
    elseif ishandle(options.ax)
        ax = options.ax;  % Use the provided axes handle if valid
    end
    
    % Adjust ColorCutMinus if it is set to -1
    if options.ColorCutMinus == -1
        options.ColorCutMinus = options.ColorCut;  % Set ColorCutMinus equal to ColorCut if it was set to -1
    end

    % Process BCCAR data by scaling and limiting values to ColorCut limits
    BCCAR = real(BCCAR) * options.scalefactor;  % Scale the BCCAR data
    maxBCplus = options.ColorCut * max(BCCAR, [], 'all');  % Maximum positive value for scaling
    maxBCMinus = options.ColorCutMinus * min(BCCAR, [], 'all');  % Maximum negative value for scaling
    maxBC = max(abs(maxBCplus), abs(maxBCMinus));  % Determine the maximum absolute value
    BCCAR(BCCAR > maxBCplus) = maxBCplus;  % Limit values above the max positive scale
    BCCAR(BCCAR < maxBCMinus) = maxBCMinus;  % Limit values below the max negative scale

    % Check if BCCAR is a matrix (not a vector), and handle accordingly
    if ~isvector(BCCAR)
        if size(Grid, 3) == 1
            % Reshape Grid into a mesh if it's 2D
            sizemesh = size(BCCAR);
            Grid1 = reshape(Grid(:, 1), sizemesh);
            Grid2 = reshape(Grid(:, 2), sizemesh);
            Grid3 = reshape(Grid(:, 3), sizemesh);
        else
            % Handle 3D grid case
            Grid1 = Grid(:, :, 1);
            Grid2 = Grid(:, :, 2);
            Grid3 = Grid(:, :, 3);
        end
        
        % Create a surface plot of BCCAR data
        h = surf(ax, Grid1, Grid2, Grid3, BCCAR, 'EdgeColor', 'none');
    else
        % If BCCAR is a vector, interpolate to create a mesh
        X = Grid(:, 1);  % X coordinates
        Y = Grid(:, 2);  % Y coordinates
        Z = Grid(:, 3);  % Z coordinates
        
        % Unique grid points for meshing
        Grid1 = unique(X);
        Grid2 = unique(Y);
        Grid3 = unique(Z);
        
        if length(Grid3) == 1
            % Handle 2D grid
            [meshX, meshY, meshZ] = meshgrid(Grid1, Grid2, Grid3);
            % Interpolate BCCAR values onto the mesh
            meshU = griddata(X, Y, BCCAR, meshX, meshY, options.method);
        else
            % Handle 3D grid
            [meshX, meshY, meshZ] = meshgrid(Grid1, Grid2, Grid3);
            % Interpolate BCCAR values onto the 3D mesh
            meshU = griddata(X, Y, Z, BCCAR, meshX, meshY, meshZ, options.method);
        end
        
        % Plot either a surface or contour based on the 'plotmode' option
        switch options.plotmode
            case 'surf'
                h = surf(ax, meshX, meshY, meshZ, meshU, 'EdgeColor', 'none');
            case 'contour'
                h = contourf(ax, meshX, meshY, meshU);  % Contour plot
        end
    end

    % Apply a custom colormap
    colormap(ax, ColorMap.redblue);

    % Reference lattice and offset for reciprocal space visualization
    Gk = (2 * pi * eye(3) / Rm).';  % Reciprocal lattice vectors
    Gknorm = norm(Gk(1, :)) / 10;  % Normalize the first reciprocal lattice vector
    offsetX = ones(2, 2) * Gknorm;  % Offset for X-direction
    offsetY = ones(2, 2) * Gknorm;  % Offset for Y-direction
    offsetY2 = ones(2, 2) * 2 * Gknorm;  % Additional offset for Y-direction

    % Get grid boundaries for plotting reference
    minx = min(Grid1, [], 'all');
    miny = min(Grid2, [], 'all');
    minz = min(Grid3, [], 'all');
    maxx = max(Grid1, [], 'all');
    maxy = max(Grid2, [], 'all');
    maxz = max(Grid3, [], 'all');
    
    % Plot boundary surfaces to indicate the limits of the grid
    h2 = surf(ax, minx - offsetX, maxy + offsetY, zeros(2), -ones(2, 2) * maxBC, 'EdgeColor', 'none');
    h3 = surf(ax, minx - offsetX, maxy + offsetY2, zeros(2), ones(2, 2) * maxBC, 'EdgeColor', 'none');
    
    % Add a colorbar for the plot
    try
        colorbar(ax, 'Ticks', [-roundn(maxBC, -2), 0, roundn(maxBC, -2)]);
    catch
        % Handle potential issues with colorbar
    end

    % Set equal axis scaling and optional shading
    axis(ax, 'equal');
    if options.shading
        shading(ax, 'interp');
    end

    % Optionally plot Brillouin Zone
    if options.BZ
        ax = BZplot(Rm, 'color', 'none', 'ax', ax, 'blackwhite', true, ...
                    'mode', options.BZmode, 'label', options.BZlabel, 'OriginPoint', options.GkOrigin);
        view(ax, 2);  % Set the view for Brillouin Zone
    end

    % Return output handles
    if nargout == 3
        varargout{1} = ax.Parent;
        varargout{2} = ax;
        varargout{3} = h;
    elseif nargout == 2
        varargout{1} = ax.Parent;
        varargout{2} = ax;
    elseif nargout == 1
        varargout{1} = ax;
    end
end
