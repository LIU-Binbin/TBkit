function varargout = quiverplot(R, A, options)
    % QUIVERPLOT Plots a 2D vector field with optional BZ plot overlay.
    %   The function displays a quiver plot of vectors based on input 
    %   positions (R) and directions (A). It can also display the Brillouin 
    %   Zone (BZ) overlay if requested.
    %
    %   Inputs:
    %       R           : Nx2 or Nx3 matrix of coordinates (position vectors).
    %       A           : Nx2 or Nx3 matrix of vectors (directions).
    %       options     : A structure containing the following optional fields:
    %           - color       : Vector color (default: 'b').
    %           - displayname : Name of the vector field (default: 'vector field').
    %           - ax          : Axes handle for the plot (optional).
    %           - scale       : Scaling factor for the vectors (default: 0.5).
    %           - ShowArrowHead : Display arrowheads (default: 'on').
    %           - Rm          : Unit cell or POSCAR data (default: read POSCAR file).
    %           - BZ          : Boolean flag to display Brillouin Zone (default: false).
    %           - BZmode      : '3D' or '2D' Brillouin Zone plot mode (default: '2D').
    %           - BZlabel     : Boolean flag to display Brillouin Zone labels (default: true).
    %
    %   Outputs:
    %       varargout   : Depending on the number of output arguments, returns the axes handles.

    arguments
        R
        A
        options.color = 'b';                      % Default color: blue
        options.displayname = 'vector field';     % Default display name
        options.ax = handle([]);                  % Optional axes handle
        options.scale = 0.5;                      % Scaling factor for vectors
        options.ShowArrowHead = 'on';             % Display arrowheads by default
        options.Rm = POSCAR_read;                 % Read POSCAR if not provided
        options.BZ logical = false;               % Show Brillouin Zone overlay (default: false)
        options.BZmode {mustBeMember(options.BZmode, {'3D', '2D'})} = '2D'; % BZ mode (2D or 3D)
        options.BZlabel logical = true;           % Show BZ labels by default
    end

    % Initialize axes
    if isempty(options.ax)
        Fig = Figs(1, 1);  % Create a new figure if no axes handle is provided
        ax = Fig.axes(1);
    else
        if ishandle(options.ax)
            ax = options.ax;  % Use the provided axes handle if valid
        else
            error('Invalid axes handle provided.');
        end
    end

    % Extract coordinates and vectors
    X = R(:, 1);
    Y = R(:, 2);
    U = A(:, 1);
    V = A(:, 2);

    % Plot the quiver (vector field)
    quiver(ax, X, Y, U, V, options.scale, 'DisplayName', options.displayname, ...
        'Color', options.color, 'LineWidth', 1, 'ShowArrowHead', options.ShowArrowHead);

    % If requested, overlay the Brillouin Zone (BZ)
    if options.BZ
        ax = TBkit_plot.BZplot(options.Rm, 'color', 'none', 'ax', ax, ...
            'blackwhite', true, 'mode', options.BZmode, 'label', options.BZlabel);
    end

    % Set the view and axis properties
    view(0, 90);  % Set the view to top-down (2D)
    axis(ax, 'equal');  % Equal scaling for x, y, and z axes

    % Return axes handles if requested
    if nargout == 2
        varargout{1} = ax.Parent;
        varargout{2} = ax;
    elseif nargout == 1
        varargout{1} = ax;
    end
end
