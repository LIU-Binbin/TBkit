function varargout = BZplot(Rm, options)
arguments
    Rm = POSCAR_read;  % Default to reading the POSCAR file if not provided
    options.ax = handle([]);  % Axis handle for plotting
    options.color = 'none';  % Default color for plot
    options.blackwhite = false;  % Whether to use black-and-white color scheme
    options.alpha = 0.2;  % Transparency for the plot
    options.mode = '3D';  % Plot mode: '3D', '2D', or '1D'
    options.Gk = false;  % Flag to use reciprocal lattice vectors (default is false)
    options.light;  % Light source for the plot (not used in current code)
    options.OriginPoint = [0 0 0];  % Origin point for the plot
    options.Rotation = eye(3);  % Rotation matrix for the plot
    options.label = true;  % Whether to label the axes and key points
    options.KPOINTS = '';  % Path to the KPOINTS file for k-point labels
    options.title = '';  % Plot title (optional)
end
import TBkit_tool_outer.*  % Import external functions from the TBkit toolset

% Extract the origin point coordinates
X0 = options.OriginPoint(1);
Y0 = options.OriginPoint(2);
Z0 = options.OriginPoint(3);

% Create a new figure and axis handle if none is provided
if isempty(options.ax)
    Fig = Figs(1, 1);
    ax = Fig.axes(1);
else
    ax = options.ax;  % Use the provided axis handle
end

% If the input Rm is an instance of the TBkit class, extract the lattice matrix
if isa(Rm, 'TBkit')
    Rm = Rm.Rm;
end

% If the lattice matrix is 2x2, extend it to 3x3 by adding a zero row and column
if isequal(size(Rm), [2, 2])
    Rm = [Rm, [0; 0]];  % Add the third dimension
end

% Set up the color and alpha (transparency) for the plot
color = options.color;
alpha = options.alpha;

% Handle the different plot modes (3D, 2D, 1D)
switch options.mode
    case '3D'
        % Calculate reciprocal lattice vectors in 3D
        if options.Gk
            Gk = Rm;  % Use the provided lattice matrix directly
        else
            Gk = (eye(3) * 2 * pi / Rm).';  % Convert to reciprocal space
        end
        if size(Gk, 1) < 3 || size(Gk, 2) < 3
            Gk(3, 3) = 1;  % Ensure 3x3 size for the lattice matrix
        end
        
        % Generate a grid of k-points for plotting
        vector = [];
        for i = linspace(-2, 2, 5)
            for j = linspace(-2, 2, 5)
                for k = linspace(-2, 2, 5)
                    vector = [vector; i, j, k];  % Generate k-points in the grid
                end
            end
        end
        
        % Rotate the reciprocal lattice vectors using the provided rotation matrix
        Gk = options.Rotation * Gk;
        vectorR = vector * Gk;  % Apply the rotation to the k-points
        
        % Find the index of the gamma point (0,0,0)
        [~, Line_000] = ismember([0, 0, 0], vectorR, 'rows');
        
        % Compute the Voronoi diagram for the k-points
        [v, c] = voronoin(vectorR);
        V = v(c{Line_000}, :);  % Extract the Voronoi vertices for the relevant region
        
        % Plot the polyhedron representing the Brillouin zone
        ax = plotPolyhedron(V, color, alpha, 'ax', ax, 'OriginPoint', options.OriginPoint);
        
        % Plot the origin (Gamma point)
        scatter3(ax, X0, Y0, Z0, 'filled', 'LineWidth', 300, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'k');
        
        % If no KPOINTS file is provided, label the Gamma point
        if strcmp(options.KPOINTS, '')
            text(ax, X0, Y0, Z0, "\Gamma", 'FontSize', 24);
        else
            % Read k-points from the KPOINTS file and label them
            [kpoints, ~, ~, kpoints_name] = KPOINTS_read(options.KPOINTS);
            [kpoints, label_unique] = unique(kpoints, 'rows');
            kpoints_name = kpoints_name(label_unique);
            
            % Replace different forms of "Gamma" with the Greek letter Γ
            kpoints_name = strrep(kpoints_name, 'GAMMA', 'Γ');
            kpoints_name = strrep(kpoints_name, 'Gamma', 'Γ');
            kpoints_name = strrep(kpoints_name, 'G', 'Γ');
            
            % Calculate the positions of the high-symmetry k-points in the reciprocal space
            HighK = kpoints * Gk + options.OriginPoint;
            
            % Plot and label the high-symmetry k-points
            for i = 1:length(HighK)
                scatter3(ax, HighK(i, 1), HighK(i, 2), HighK(i, 3), 'filled', 'LineWidth', 300);
                text(ax, HighK(i, 1), HighK(i, 2), HighK(i, 3), kpoints_name(i), 'FontSize', 24);
            end
        end
        
        % Label the axes
        xlabel('k_x');
        ylabel('k_y');
        zlabel('k_z');
        
        % If labels are enabled, plot the reciprocal lattice vectors
        if options.label
            if options.blackwhite
                % Plot black-and-white quivers for the reciprocal lattice vectors
                quiver3(ax, X0, Y0, Z0, Gk(1, 1), Gk(1, 2), Gk(1, 3), 'Color', 'k', 'DisplayName', 'k_1', 'AutoScale', 'off');
                quiver3(ax, X0, Y0, Z0, Gk(2, 1), Gk(2, 2), Gk(2, 3), 'Color', 'k', 'DisplayName', 'k_2', 'AutoScale', 'off');
                quiver3(ax, X0, Y0, Z0, Gk(3, 1), Gk(3, 2), Gk(3, 3), 'Color', 'k', 'DisplayName', 'k_3', 'AutoScale', 'off');
            else
                % Plot color-coded quivers for the reciprocal lattice vectors
                quiver3(ax, X0, Y0, Z0, Gk(1, 1), Gk(1, 2), Gk(1, 3), 'Color', 'r', 'DisplayName', 'k_1', 'AutoScale', 'off');
                quiver3(ax, X0, Y0, Z0, Gk(2, 1), Gk(2, 2), Gk(2, 3), 'Color', 'g', 'DisplayName', 'k_2', 'AutoScale', 'off');
                quiver3(ax, X0, Y0, Z0, Gk(3, 1), Gk(3, 2), Gk(3, 3), 'Color', 'b', 'DisplayName', 'k_3', 'AutoScale', 'off');
            end
            
            % Shift the lattice vectors and label them
            Gk = Gk + options.OriginPoint;
            text(ax, Gk(1, 1), Gk(1, 2), Gk(1, 3), "k_1", 'FontSize', 24);
            text(ax, Gk(2, 1), Gk(2, 2), Gk(2, 3), "k_2", 'FontSize', 24);
            text(ax, Gk(3, 1), Gk(3, 2), Gk(3, 3), "k_3", 'FontSize', 24);
        end
        
    case '2D'
        % Handle the case for a 2D Brillouin zone
        % Similar to the 3D case, but now reduced to two dimensions
        if size(Rm, 1) == 3
            Rm = Rm(1:2, :);  % Use only the 2D portion of the matrix
        end
        if options.Gk
            Gk = Rm;  % Use the provided lattice matrix directly
        else
            Gk = (eye(3) * 2 * pi / Rm).';  % Convert to reciprocal space
        end
        
        % Apply the rotation to the reciprocal lattice vectors
        Gk = Gk * options.Rotation;
        
        % Convert the reciprocal lattice vectors to a 2D basis
        Gk_x = Gk(1, :) / norm(Gk(1, :));
        Gk_z = cross(Gk_x, Gk(2, :));
        Gk_z = Gk_z / norm(Gk_z);
        Gk_y = cross(Gk_z, Gk_x);
        Gk_2d_Pmat = [Gk_x; Gk_y];
        
        % Calculate the positions of the k-points in the 2D reciprocal space
        Gk_2d = Gk * Gk_2d_Pmat.';
        vector = [];
        for i = linspace(-3, 3, 7)
            for j = linspace(-3, 3, 7)
                vector = [vector; i, j];
            end
        end
        
        % Apply the 2D transformation to the vector list
        vectorR = vector * Gk_2d;
        [~, Line_000] = ismember([0, 0], vectorR, 'rows');
        [v, c] = voronoin(vectorR);
        V = v(c{Line_000}, :);
        
        % Plot the 2D Brillouin zone
        ax = plotPolyhedron(V, color, alpha, 'ax', ax, 'OriginPoint', options.OriginPoint, 'Orientation', Gk_2d_Pmat);
        
        % Plot the origin (Gamma point) and label it
        scatter3(ax, X0, Y0, Z0, 'filled', 'LineWidth', 300);
        text(ax, X0, Y0, Z0, "\Gamma", 'FontSize', 24);
        
        % Label the axes
        xlabel('k_x');
        ylabel('k_y');
        zlabel('k_z');
        
        % If labels are enabled, plot the reciprocal lattice vectors
        if options.label
            quiver3(ax, X0, Y0, Z0, Gk(1, 1), Gk(1, 2), Gk(1, 3), 'Color', 'r', 'DisplayName', 'k_1', 'AutoScale', 'off');
            quiver3(ax, X0, Y0, Z0, Gk(2, 1), Gk(2, 2), Gk(2, 3), 'Color', 'g', 'DisplayName', 'k_2', 'AutoScale', 'off');
            view(ax, 0, 90);  % Top-down view for 2D plot
        end
        
    case '1D'
        % Handle the case for a 1D Brillouin zone
        % In this case, only 1D lattice vectors would be plotted
        % This case can be extended as needed.
end

% Return the plot handles if requested
if nargout == 2
    varargout{1} = ax.Parent;
    varargout{2} = ax;
end
if nargout == 1
    varargout{1} = ax;
end
end
