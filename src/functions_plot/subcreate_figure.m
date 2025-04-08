function [fig, ax] = subcreate_figure(m, n, position, options)
    % SUBCREATE_FIGURE Creates a figure with subplots based on the given arguments.
    %   This function creates a figure with `m` rows and `n` columns of subplots.
    %
    %   Inputs:
    %       m        : Number of rows for the subplot grid (default: 2).
    %       n        : Number of columns for the subplot grid (default: 1).
    %       position : Position of the figure (default: auto based on `m` and `n`).
    %       options  : Struct with optional fields:
    %           - papersize  : Scaling factor for the figure (default: 1).
    %           - fontname   : Font for axes (default: 'Helvetica').
    %           - color      : Background color of the figure (default: 'white').
    %           - fig        : Existing figure handle (default: create new).
    %           - ax         : Existing axes handles (default: create new).
    %           - FontSize   : Font size for axes (default: 12).
    %
    %   Outputs:
    %       fig      : Figure handle.
    %       ax       : Array of axes handles for the subplots.

    arguments
        m double {mustBeInteger} = 2;             % Number of rows for subplots (default 2)
        n double {mustBeInteger} = 1;             % Number of columns for subplots (default 1)
        position double = [];                     % Position of the figure (optional)
        options.papersize double = 1;             % Scaling factor for the figure (default 1)
        options.fontname string = "Helvetica";    % Font name for axes (default 'Helvetica')
        options.color string = 'white';           % Background color of the figure (default 'white')
        options.fig handle = [];                  % Existing figure handle (optional)
        options.ax handle = [];                   % Existing axes handle (optional)
        options.FontSize double = 12;             % Font size for axes (default 12)
    end

    % Determine figure position if not provided
    if isempty(position)
        maxmn = max(m, n);
        if n == maxmn
            height = 0.8 * m / n;
            width = 0.6;
        elseif m == maxmn
            height = 0.8;
            width = 0.6 * n / m;
        end
        position = [0.1 0.1 width height];
    end
    position([3, 4]) = position([3, 4]) * options.papersize;

    % Create the figure or use the provided one
    if isempty(options.fig)
        fig = figure('Unit', 'normalized', 'Position', position, 'Color', options.color);
    else
        fig = options.fig;
    end

    % Create the axes or use the provided ones
    if isempty(options.ax)
        for i = 1:(m * n)
            ax(i) = subplot(m, n, i, 'Parent', fig, 'LineWidth', 1.5, ...
                            'FontSize', options.FontSize, 'FontName', options.fontname);
        end
    else
        ax = options.ax;
    end

    % Set box and hold for all axes
    for i = 1:length(ax)
        box(ax(i), 'on');
        hold(ax(i), 'all');
    end
end
