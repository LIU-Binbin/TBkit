function varargout = WilsonLoopPlot(BFCAR, klist_l, options)
    % WILSONLOOPPLOT Plots the Wilson loop data from BFCAR data and given k-points.
    % This function visualizes the Berry phase (or Wannier center) across a given 
    % k-space path. It allows interactive customization of the plot.
    %
    % Inputs:
    %   BFCAR          - Berry phase (or Wannier center) data (matrix form).
    %   klist_l        - List of k-points for the plot.
    %   options        - Optional parameters to control the plot appearance and behavior.
    %       ycut      - Y-axis range for the plot [default: [-pi, pi]].
    %       shift     - Shift the data by this value [default: 0].
    %       Color     - RGB values for plot color [default: random color].
    %       LineSpec  - Line specification for plotting (e.g., 'o', '-', etc.) [default: 'o'].
    %       title     - Title of the plot [default: 'WilsonLoop'].
    %       xlabel    - X-axis label [default: 'k_{envolution}'].
    %       ylabel    - Y-axis label [default: 'WannierCenter'; 'BerryPhase1D'].
    %       ax        - Axes handle to plot on (if provided).
    %
    % Outputs:
    %   varargout     - If requested, returns the axis or the parent of the axis.
    
    % Parse and set options with default values if necessary.
    arguments
        BFCAR;
        klist_l;
        options.ycut = [-pi, pi];             % Default Y-axis cut range.
        options.shift = 0;                     % Default shift for data.
        options.Color = [rand rand rand];      % Random color by default.
        options.LineSpec = 'o';                % Default line specification.
        options.title = 'WilsonLoop';          % Default plot title.
        options.xlabel = 'k_{envolution}';    % Default X-axis label.
        options.ylabel = ["WannierCenter"; "(BerryPhase1D)"];  % Default Y-axis label.
        options.ax = handle([]);               % Default is empty, meaning new axes will be created.
    end
    
    % Adjust the BFCAR values by shifting negative values by 2*pi (wrap-around for periodic data).
    if options.ycut(1) >= 0
        BFCAR(BFCAR(:,:) < 0) = 2 * pi + BFCAR(BFCAR(:,:) < 0);
    end
    
    % Create new axes if none are provided in options.
    if isempty(options.ax)
        [~, ax] = Figs;  % Create a new figure if no axes are specified.
    else
        ax = options.ax;  % Use provided axes handle.
    end
    
    % Ensure klist_l has two values, if not, extend it.
    if length(klist_l) == 1
        BFCAR = [BFCAR, BFCAR];  % Duplicate BFCAR data for plotting.
        klist_l = [0, klist_l];   % Extend k-points for plotting.
    end
    
    % Plot the Wilson loop (Berry phase/Wannier center) with the given options.
    ax = bandplot(BFCAR + options.shift, options.ycut + options.shift, klist_l, [], [], ...
        'LineSpec', options.LineSpec, ...
        'Color', options.Color, ...
        'xlabel', options.xlabel, ...
        'ylabel', options.ylabel, ...
        'ax', ax, ...
        'title', options.title);
    
    % Customize the x and y ticks with specific values for periodicity visualization.
    set(ax, 'xtick', -2*pi:pi/2:2*pi, ...
            'xticklabel', ["-2\pi", "-3\pi/2", "-\pi", "-\pi/2", ...
                          "0", "\pi/2", "\pi", "3\pi/2", "2\pi"]);
    set(ax, 'ytick', -2*pi:pi/2:2*pi, ...
            'yticklabel', ["-2\pi", "-3\pi/2", "-\pi", "-\pi/2", ...
                          "0", "\pi/2", "\pi", "3\pi/2", "2\pi"]);
    
    % Return output depending on the number of requested outputs.
    if nargout == 2
        varargout{1} = ax.Parent;  % Return the parent figure of the axes.
        varargout{2} = ax;         % Return the axes handle.
    end
    
    if nargout == 1
        varargout{1} = ax;         % Return the axes handle only.
    end

end
