function varargout = WccPlot(WccCAR, klist_l, options)
    % WCCPLOT Plots the Wcc data (e.g., Wannier Charge Center) as a function of k-points.
    % This function visualizes Wcc data across a given k-space path. It allows customization
    % of the plot appearance, such as color, axis labels, and title.
    %
    % Inputs:
    %   WccCAR        - Wcc data (matrix form, for each k-point).
    %   klist_l       - List of k-points corresponding to the Wcc data.
    %   options       - Optional parameters to control the plot appearance and behavior.
    %       ycut      - Y-axis range for the plot [default: [0, 1]].
    %       shift     - Shift the data by this value [default: 0].
    %       Color     - RGB values for plot color [default: random color].
    %       LineSpec  - Line specification for plotting (e.g., 'o', '-', etc.) [default: 'o'].
    %       title     - Title of the plot [default: ''].
    %       xlabel    - X-axis label [default: 'k_{envolution}'].
    %       ylabel    - Y-axis label [default: 'Wcc'].
    %       ax        - Axes handle to plot on (if provided).
    %
    % Outputs:
    %   varargout     - If requested, returns the axis or the parent of the axis.
    
    % Parse and set options with default values if necessary.
    arguments
        WccCAR;
        klist_l;
        options.ycut = [0, 1];              % Default Y-axis cut range.
        options.shift = 0;                   % Default shift for data.
        options.Color = [rand rand rand];    % Random color by default.
        options.LineSpec = 'o';              % Default line specification.
        options.title = '';                  % Default plot title (empty).
        options.xlabel = 'k_{envolution}';  % Default X-axis label.
        options.ylabel = ["Wcc"];            % Default Y-axis label.
        options.ax = handle([]);             % Default is empty, meaning new axes will be created.
    end
    
    % Adjust the WccCAR values for proper periodicity handling.
    % If negative values are present and ycut(1) >= 0, shift them.
    if options.ycut(1) >= 0
        WccCAR(WccCAR(:,:) < 0) = 1 + WccCAR(WccCAR(:,:) < 0);  % Adjust negative values for periodicity.
    end
    
    % Create new figure and axes if none are provided in options.
    if isempty(options.ax)
        Fig = create_figure('Position', [0.2, 0.2, 0.6, 0.6]);  % Create a figure with specific position.
        ax = Fig.axes(1);  % Get the axes from the figure.
    else
        ax = options.ax;  % Use the provided axes handle if available.
    end
    
    % Ensure klist_l has two values. If only one is provided, duplicate data for plotting.
    if length(klist_l) == 1
        WccCAR = [WccCAR, WccCAR];  % Duplicate Wcc data for k-points.
        klist_l = [0, klist_l];     % Extend k-points for the plot.
    end
    
    % Plot the Wcc data using the bandplot function with the given options.
    ax = bandplot(WccCAR + options.shift, options.ycut + options.shift, klist_l, [], [], ...
        'LineSpec', options.LineSpec, ...
        'Color', options.Color, ...
        'xlabel', options.xlabel, ...
        'ylabel', options.ylabel, ...
        'ax', ax, ...
        'title', options.title);
    
    % Customize the x and y ticks for periodicity and visualization.
    set(ax, 'xtick', -2*pi:pi/2:2*pi, ...
            'xticklabel', ["-2\pi", "-3\pi/2", "-\pi", "-\pi/2", ...
                          "0", "\pi/2", "\pi", "3\pi/2", "2\pi"]);
    set(ax, 'ytick', -1:0.25:1, ...
            'yticklabel', ["-1", "-0.75", "-0.5", "-0.25", ...
                          "0", "0.25", "0.5", "0.75", "1"]);
    
    % Return output depending on the number of requested outputs.
    if nargout == 2
        varargout{1} = ax.Parent;  % Return the parent figure of the axes.
        varargout{2} = ax;         % Return the axes handle.
    end
    
    if nargout == 1
        varargout{1} = ax;         % Return the axes handle only.
    end

end
