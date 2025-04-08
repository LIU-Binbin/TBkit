function varargout = set_reference(kpoints_l, kpoints_name, X_cut, Y_cut, options)
% set_reference: A function to customize the plotting of band structures.
%
% This function allows setting up the x-axis for k-points with their labels, 
% adjusting plot limits, and customizing font and axis labels.
%
% Syntax:
%   [fig, ax] = set_reference(kpoints_l, kpoints_name, X_cut, Y_cut, options)
%   [fig, ax] = set_reference(kpoints_l, kpoints_name, X_cut, Y_cut)
%   [fig, ax] = set_reference(kpoints_l, kpoints_name)
%
% Input:
%   kpoints_l   - Numeric array representing the positions of k-points on the x-axis.
%   kpoints_name- Cell array of strings containing the names or labels of k-points.
%   X_cut       - 2-element vector specifying x-axis limits [X_min, X_max].
%   Y_cut       - 2-element vector specifying y-axis limits [Y_min, Y_max].
%   options     - Structure with optional fields for customizing the plot:
%     .mode      - Plot mode ('band', default).
%     .ax        - Axis handle to plot on (default: current axis).
%     .FontName  - Font name for labels (default: 'Helvetica').
%     .FontSize  - Font size for labels (default: 24).
%     .xlabel    - Custom x-axis label (default: '').
%     .ylabel    - Custom y-axis label (default: 'E(eV)').
%
% Output:
%   ax          - The axis handle used for the plot.
%   fig         - The figure handle containing the plot.
%
% Example:
%   kpoints_l = [0, 0.5, 1];
%   kpoints_name = {'\Gamma', 'X', 'M'};
%   X_cut = [0, 1];
%   Y_cut = [-5, 5];
%   options = struct('mode', 'band', 'FontName', 'Arial', 'FontSize', 16, 'xlabel', 'k-points', 'ylabel', 'Energy (eV)');
%   [fig, ax] = set_reference(kpoints_l, kpoints_name, X_cut, Y_cut, options);
%
% Change Log:
%   Document Date: 2020/12/03
%   Creation Date: 2020/12/03
%   Last updated: 2020/12/03
%
% Copyright:
%   parkman
%   <parkman@buaa.edu.cn>

%--------  Initialization --------
arguments 
    kpoints_l = [];
    kpoints_name = [];
    X_cut = [];
    Y_cut = [];
    options.mode = 'band'; % Default mode is 'band'
    options.ax = handle([]); % Default axis is current axis
    options.FontName = 'Helvetica';
    options.FontSize = 24;
    options.xlabel = '';
    options.ylabel = 'E(eV)';
end

% Set font name
FontName = options.FontName;

%--------  Axis handle setup --------
if isempty(options.ax)
    ax = gca; % Use current axis if not provided
else
    ax = options.ax;
end

%--------  Main Plotting Logic --------
if strcmp(options.mode, 'band')
    % Check and set x and y limits
    try
        set(ax, 'XLim', X_cut);
    catch
        warning('Invalid X_cut. Check KPOINTS or bandplot settings.');
    end
    set(ax, 'YLim', Y_cut);
    
    % Set font, label, and ticks
    set(ax, 'FontName', options.FontName, 'FontSize', options.FontSize, 'LineWidth', 1, ...
        'XTick', kpoints_l, 'XTickLabel', kpoints_name);
    ylabel(ax, options.ylabel, 'FontName', FontName);
    xlabel(ax, options.xlabel, 'FontName', FontName);
    
    % Plot vertical lines for each k-point
    for i = 1:(length(kpoints_name) - 2)
        X = [kpoints_l(i+1), kpoints_l(i+1)];
        plot(ax, X, Y_cut, 'LineWidth', 0.1, 'Color', [0, 0, 0], 'DisplayName', 'K-path', ...
            'HandleVisibility', 'off');
    end
end

% Hold the plot for further modifications
hold(ax, 'on');
axis(ax, 'square'); % Square the axes to maintain equal scaling

%--------  Return Output --------
if nargout == 2
    varargout{1} = ax.Parent; % Return figure handle
    varargout{2} = ax;        % Return axis handle
elseif nargout == 1
    varargout{1} = ax;        % Return axis handle
end

end
