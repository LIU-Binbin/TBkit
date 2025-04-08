function varargout = braidplot1D(EIGENCAR, Ecut, klist_l, kpoints_l, kpoints_name, options)
arguments
    EIGENCAR = EIGENVAL_read();  % Default to reading EIGENVAL if not provided
    Ecut double = [-3, 3; -3, 3];  % Energy range for plotting
    klist_l double = [];  % k-point list
    kpoints_l double = [];  % k-point list in numerical form
    kpoints_name string = [];  % k-point names
    options.ax = handle([]);  % Axis handle for the plot (default: new figure)
    options.FontName = 'Helvetica';  % Font for labels
    options.FontSize = 24;  % Font size for labels
    options.KPOINTS = 'KPOINTS';  % Default KPOINTS file
    options.POSCAR = 'POSCAR';  % Default POSCAR file
    options.Color = @parula;  % Default colormap (parula)
    options.title = '';  % Title of the plot (empty string if none)
    options.klist_l = [];  % klist for plotting
    options.kpoints_l = [];  % k-points labels
    options.kpoints_name = [];  % k-point names for labeling
    options.xlabel = '';  % x-axis label
    options.ylabel = 'E (eV)';  % y-axis label
    options.LineSpec = '-';  % Line style for plot
    options.LineWidth = 3;  % Line width for plot
    options.MarkerSize = 3;  % Marker size
    options.MarkerEdgeColor = 'none';  % Marker edge color
    options.MarkerFaceColor = 'none';  % Marker face color
    options.set_reference = true;  % Set reference (grid and labels) on the plot
    options.set_plane = true;  % Set plane markers on the plot
    options.Units = 'pixels';  % Units for plot sizing
    options.Position = [];  % Position for the figure
    options.flip = true;  % Flip the colormap
    options.Mirror = -1;  % Mirror the plot in the z-direction
end

%-------- Initialize the function --------
import TBkit_tool.*  % Import tool functions from TBkit
nOutputs = nargout;  % Get number of output arguments
varargout = cell(1, nOutputs);  % Initialize output cell array

%-------- Handle missing arguments --------
if isempty(klist_l) 
    Rm = POSCAR_read(options.POSCAR);  % Read the lattice matrix from POSCAR
    [~, klist_l, ~, kpoints_l, kpoints_name] = kpathgen3D(Rm, options.KPOINTS);  % Generate k-point path
end

% Get number of bands (rows in EIGENCAR)
Nbands = size(EIGENCAR, 1);

% Set up colormap based on options.Color
if isstring(options.Color)
    if isrow(options.Color)
        colorL = options.Color.';  % Transpose the color string if needed
    end
elseif isnumeric(options.Color)
    x = linspace(0, 1, size(options.Color, 1));  % Linearly spaced x values for color interpolation
    xq = linspace(0, 1, Nbands);  % Interpolated values for Nbands
    colorL = [
        interp1(x, options.Color(:,1), xq).',  % Interpolated red channel
        interp1(x, options.Color(:,2), xq).',  % Interpolated green channel
        interp1(x, options.Color(:,3), xq).',  % Interpolated blue channel
    ];
elseif isa(options.Color, 'function_handle')
    colorL = options.Color(Nbands);  % Apply function handle if passed as colormap
else
    for i = 1:Nbands
        colorL(i,:) = [rand, rand, rand];  % Random colors if none specified
    end
end

% Flip the colormap if the option is set
if options.flip
    colorL = flip(colorL, 1);
end

% Set the marker edge and face colors based on options
if isempty(options.MarkerEdgeColor)
    MarkerEdgeColor = colorL;
else
    MarkerEdgeColor = options.MarkerEdgeColor;
end

if isempty(options.MarkerFaceColor)
    MarkerFaceColor = colorL;
else
    MarkerFaceColor = options.MarkerFaceColor;
end

% Set the plot title based on the current directory or user-specified title
if strcmp(options.title, '')
    [~, titlestring] = fileparts(pwd);  % Use the current folder name as title
    titlestring = strrep(titlestring, '_', '-');  % Replace underscores with hyphens
else
    titlestring = options.title;  % Use user-specified title
end

% Create a new figure if no axis is passed, or use the provided axis handle
if isempty(options.ax)
    Fig = Figs(1, 1);  % Create a new figure
    ax = Fig.axes(1);  % Get the axis handle
else
    if isvalid(options.ax)
        ax = options.ax;  % Use the provided axis if valid
    else
        Fig = Figs(1, 1);
        ax = Fig.axes(1);
    end
end

%-------- Handle k-point list and names --------
if isempty(options.klist_l)
    klist = klist_l;  % Use the provided klist
else
    klist = options.klist_l;  % Use the option klist if set
end

if isempty(options.kpoints_l)
    % kpoints_l = kpoints_l;  % Not needed as it's passed through options
else
    kpoints_l = options.kpoints_l;
end

if isempty(options.kpoints_name)
    % kpoints_name = kpoints_name;  % Not needed as it's passed through options
else
    kpoints_name = options.kpoints_name;
end

% Determine the x-axis limits based on the k-point list
xmin = klist(1);
xmax = klist(end);
Xcut = [xmin, xmax];

%-------- Set reference plane and plot settings --------
if options.set_reference
    % Reference settings for grid and labels
    Ycut = Ecut(1, :);  % Set y-axis energy range
    Zcut = Ecut(2, :);  % Set z-axis energy range
    axis(ax, 'square');  % Make the axis square
    ylabel(ax, 'Re(E(z))');  % Label for the real part of energy
    zlabel(ax, 'Im(E(z))');  % Label for the imaginary part of energy
    xlabel(ax, 'kpath');  % Label for the k-path
    
    % Set k-point ticks and labels
    xticks(ax, kpoints_l);
    xticklabels(ax, kpoints_name);
    xlim(ax, Xcut);
    ylim(ax, Ycut);
    zlim(ax, Zcut);
    
    % Draw reference planes at each k-point if options.set_plane is true
    if options.set_plane
        nkl = length(kpoints_l);  % Number of k-points
        for i = 1:nkl
            V{i} = [
                kpoints_l(i), Ycut(1), Zcut(1);
                kpoints_l(i), Ycut(1), Zcut(2);
                kpoints_l(i), Ycut(2), Zcut(2);
                kpoints_l(i), Ycut(2), Zcut(1);
            ];
            patch(ax, 'Faces', [1, 2, 3, 4], 'Vertices', V{i}, ...
                'EdgeColor', 'none', 'FaceColor', [102 102 102] / 255, 'FaceAlpha', 0.1, ...
                'DisplayName', ['T_', num2str(i)]);  % Draw patch for each plane
        end
    end
end

%-------- Plot the bands --------
RealEIGENCAR = real(EIGENCAR);  % Real part of EIGENCAR
ImagEIGENCAR = imag(EIGENCAR);  % Imaginary part of EIGENCAR
nKinEIGENCAR = size(EIGENCAR, 2);  % Number of k-points in EIGENCAR

% Plot the energy bands for each eigenstate
for Ei = 1:Nbands
    plot3(ax, klist(1:nKinEIGENCAR), RealEIGENCAR(Ei, :), options.Mirror * ImagEIGENCAR(Ei, :), options.LineSpec, ...
        'LineWidth', options.LineWidth, ...
        'Color', colorL(Ei, :), ...
        'MarkerSize', options.MarkerSize, ...
        'MarkerEdgeColor', MarkerEdgeColor, ...
        'MarkerFaceColor', MarkerFaceColor, ...
        'DisplayName', ['E_', num2str(Ei)]);  % Plot each energy band
end

%-------- Add title --------
title(ax, titlestring);  % Add title to the plot

%-------- Save figure if needed --------
if length(varargout) > 2
    % Save figure and associated data
    plot_data_coll.EIGENCAR = EIGENCAR;
    plot_data_coll.Ecut = Ecut;
    plot_data_coll.titlestring = titlestring;
    plot_data_coll.color = colorL;
    plot_data_coll.klist_l = klist_l;
    plot_data_coll.kpoints_l = kpoints_l;
    plot_data_coll.kpoints_name = kpoints_name;
    plot_data_coll.fontname = options.fontname;
    [~, ax] = save_figure(ax.Parent.Parent, titlestring + ".eps", ax);
    varargout{3} = plot_data_coll;
end

%-------- Return output --------
if nargout == 2
    varargout{1} = ax.Parent;
    varargout{2} = ax;
end

if nargout == 1
    varargout{1} = ax;
end
end
