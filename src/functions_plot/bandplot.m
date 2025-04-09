function varargout = bandplot(EIGENCAR, Ecut, klist_l, kpoints_l, kpoints_name, options)
% BANDPLOT Plot band structure from eigenvalues.
% This function visualizes the band structure by plotting the eigenvalues
% given in EIGENCAR against kpoints, with configurable options for
% appearance and layout.
%
% Inputs:
%   EIGENCAR     - Eigenvalue data (NxM matrix or cell array)
%   Ecut         - Energy cut range (2-element vector)
%   klist_l      - List of k-points (optional)
%   kpoints_l    - K-point labels (optional)
%   kpoints_name - Names of k-points (optional)
%   options      - Structure specifying plot options including:
%     ax             - Axes handle (optional)
%     FontName       - Font for labels (default 'Helvetica')
%     FontSize       - Font size for text (default 24)
%     KPOINTS        - KPOINTS file name (default 'KPOINTS')
%     POSCAR         - POSCAR file name (default 'POSCAR')
%     Color          - Color scheme for bands (default @jet)
%     title          - Plot title (default '')
%     xlabel         - X-axis label (default '')
%     ylabel         - Y-axis label (default 'E (eV)')
%     LineSpec       - Line style for plots (default '-')
%     LineWidth      - Line width for plots (default 1)
%     MarkerSize     - Size of markers (default 3)
%     MarkerEdgeColor - Color for marker edges (default depends on color)
%     MarkerFaceColor - Color for marker faces (default depends on color)
%     set_reference  - Boolean to set reference lines (default true)
%     Units          - Units for the figure (default 'pixels')
%     Position       - Position for the figure (default [])
%     legends        - Legends for the plots (default [])
%
% Outputs:
%   varargout      - Optional output can include axes handles
%                    and additional plot data.

arguments
    EIGENCAR = EIGENVAL_read(); % Default eigenvalue data
    Ecut double = [-3, 3]; % Default energy cut range
    klist_l double = []; % Default k-point list
    kpoints_l double = []; % Default k-points
    kpoints_name string = []; % Default k-point names
    options.ax = handle([]); % Axes handle option
    options.FontName = 'Helvetica'; % Default font
    options.FontSize = 24; % Default font size
    options.KPOINTS = 'KPOINTS'; % Default KPOINTS file
    options.POSCAR = 'POSCAR'; % Default POSCAR file
    options.Color = @jet; % Default color scheme
    options.title = ''; % Default title
    options.klist_l = []; % Default K-point list
    options.kpoints_l = []; % Default K-point positions
    options.kpoints_name = []; % Default K-point names
    options.xlabel = ''; % Default x-label
    options.ylabel = 'E (eV)'; % Default y-label
    options.LineSpec = '-'; % Default line specification
    options.LineWidth = 1; % Default line width
    options.MarkerSize = 3; % Default marker size
    options.MarkerEdgeColor = []; % Default marker edge color
    options.MarkerFaceColor = []; % Default marker face color
    options.set_reference = true; % Default reference line setting
    options.Units = 'pixels'; % Default units
    options.Position = []; % Default position
    options.legends = []; % Default legends
end

%-------- Initialize Variables --------
import TBkit_tool.*; % Import required tools
nOutputs = nargout; % Get number of outputs
varargout = cell(1, nOutputs); % Initialize output cell array

%-------- Handle K-point List and Read Files --------
if isempty(klist_l)
    Rm = POSCAR_read(options.POSCAR); % Read structure from POSCAR
    [~, klist_l, ~, kpoints_l, kpoints_name] = kpathgen3D(Rm, options.KPOINTS); % Generate k-point paths
end
%-------- Determine Number of Bands and Color Options --------
if iscell(EIGENCAR)
    Nbands = size(EIGENCAR{1}, 1); % Number of bands from cell array
    color = options.Color; % Use specified color scheme
else
    Nbands = size(EIGENCAR, 1); % Number of bands from matrix
    if ischar(options.Color) || isnumeric(options.Color)
        color = options.Color; % Use specified color if it's a char or numeric
    else
        color = [rand, rand, rand]; % Random color if no valid color specified
    end
end

% Determine marker colors
if isempty(options.MarkerEdgeColor)
    MarkerEdgeColor = color; % Default edge color if not specified
else
    MarkerEdgeColor = options.MarkerEdgeColor; % Use specified edge color
end

if isempty(options.MarkerFaceColor)
    MarkerFaceColor = color; % Default face color if not specified
else
    MarkerFaceColor = options.MarkerFaceColor; % Use specified face color
end

% Set plot title
if strcmp(options.title, '')
    [~, titlestring] = fileparts(pwd); % Use current directory name as title
    titlestring = strrep(titlestring, '_', '-'); % Replace underscores with dashes
else
    titlestring = options.title; % Use specified title
end

% Prepare axes for plotting
if isempty(options.ax)
    [Fig,ax]= Figs(1, 1); % Create new figure
else
    if isvalid(options.ax)
        ax = options.ax; % Use provided axis handle
    else
        [Fig,ax]= Figs(1, 1); % Create new figure 
    end
end

%-------- Determine K-point List and Names --------
if isempty(options.klist_l)
    klist = klist_l; % Use provided k-point list
else
    klist = options.klist_l; % Override with specified list if available
end

if isempty(options.kpoints_l)
    % Preserve original kpoints if options.kpoints_l is empty
else
    kpoints_l = options.kpoints_l; % Use specified kpoints
end

if isempty(options.kpoints_name)
    % Preserve original names if options.kpoints_name is empty
else
    kpoints_name = options.kpoints_name; % Use specified kpoint names
end

% Set X limits for the plot
xmin = klist(1);
xmax = klist(end);
Xcut = [xmin, xmax]; % Define limits for X-axis

if options.set_reference
    %-------- Set Reference Lines on the Axes --------
    [ax] = set_reference(kpoints_l, kpoints_name, Xcut, Ecut,...
        'ax', ax,...
        'fontname', options.FontName,...
        'FontSize', options.FontSize,...
        'xlabel', options.xlabel,...
        'ylabel', options.ylabel ...
        );
end

%-------- Plot Band Structure --------
if iscell(EIGENCAR)
    Npara = length(EIGENCAR); % Number of datasets in cell array

    if isa(color, 'function_handle')
        Colormap = color(Npara); % Get colormap for each dataset
    else
        Colormap = color; % Use specified color
    end

    % Plot each dataset separately
    for j = 1:Npara
        LineWidth = options.LineWidth; % Default line width

        if length(options.LineWidth) == Npara
            LineWidth = options.LineWidth(j); % Override for specific dataset
        end

        Nbands = size(EIGENCAR{j}, 1); % Number of bands for the current dataset
        for Ei = 1:Nbands
            line = plot(ax, klist, EIGENCAR{j}(Ei, :), options.LineSpec,...
                'LineWidth', LineWidth,...
                'Color', Colormap(j, :),...
                'DisplayName', num2str(j) + "_" + num2str(Ei));

            if Ei ~= 1
                line.HandleVisibility = 'off'; % Hide handles for multiple bands
            end

            hold(ax, 'on'); % Retain current plot
        end
    end
    legend(ax, options.legends); % Add legend if specified

else
    % Single set of eigenvalues
    if isa(color, 'function_handle')
        color = color(1); % Use single color for plot
    end

    % Plot each band
    for Ei = 1:Nbands
        plot(ax, klist, EIGENCAR(Ei, :), options.LineSpec,...
            'LineWidth', options.LineWidth,...
            'Color', color,...
            'MarkerSize', options.MarkerSize,...
            'MarkerEdgeColor', MarkerEdgeColor,...
            'MarkerFaceColor', MarkerFaceColor,...
            'DisplayName', num2str(Ei));
    end
end
%-------- Set Title of the Plot --------
if isa(titlestring, 'cell')
    titlestring = cell2mat(titlestring); % Convert cell to char array if needed
end
title(ax, titlestring,'FontWeight','normal'); % Set the title of the plot

%-------- Save Additional Output Data if Requested --------
if length(varargout) > 2
    % Save plot data for further analysis or export
    plot_data_coll.EIGENCAR = EIGENCAR; % Store eigenvalue data
    plot_data_coll.Ecut = Ecut; % Store energy cut information
    plot_data_coll.titlestring = titlestring; % Store title
    plot_data_coll.color = color; % Store color information
    plot_data_coll.klist_l = klist_l; % Store klist
    plot_data_coll.kpoints_l = kpoints_l; % Store kpoints
    plot_data_coll.kpoints_name = kpoints_name; % Store kpoints names
    plot_data_coll.fontname = options.FontName; % Store font info

    % Generate plot and save as EPS file
    [~, ax] = save_figure(ax.Parent, titlestring + ".eps", ax); % Save figure in EPS format
    varargout{3} = plot_data_coll; % Return additional plot data
end

%-------- Return Outputs --------
if nargout == 2
    varargout{1} = ax.Parent; % Return parent figure handle if requested
    varargout{2} = ax; % Return axis handle if requested
elseif nargout == 1
    varargout{1} = ax; % Return axis handle if only one output is requested
end
end
