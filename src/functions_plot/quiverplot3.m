function varargout = quiverplot3(R, A, color, displayname, ax, scale)
    % QUIVERPLOT3 Generates a 3D quiver plot for vector fields.
    %   Inputs:
    %       R: Nx3 matrix of the starting points (x, y, z) for the vectors.
    %       A: Nx3 matrix of vector components (U, V, W) at each point.
    %       color: Color for the quiver plot (default is 'b').
    %       displayname: Name for the display legend (default is 'vector field').
    %       ax: Axis handle where the plot will be drawn (optional).
    %       scale: Scaling factor for the vectors (default is 0.5).
    %   Outputs:
    %       If nargout >= 1, returns axis handles.

    % Set default argument values using the "arguments" block
    arguments
        R (:, 3) double  % R should be a Nx3 matrix (vector positions)
        A (:, 3) double  % A should be a Nx3 matrix (vector components)
        color (1, :) char = 'b'  % Default color is 'b' (blue)
        displayname (1, :) char = 'vector field'  % Default displayname
        ax = gca  % Default axis is the current axis (gca)
        scale double = 0.5  % Default scale factor
    end
    
    % Extract components of the positions and vectors
    X = R(:, 1);  % X coordinates of the vector origins
    Y = R(:, 2);  % Y coordinates of the vector origins
    Z = R(:, 3);  % Z coordinates of the vector origins
    U = A(:, 1);  % X component of the vector
    V = A(:, 2);  % Y component of the vector
    W = A(:, 3);  % Z component of the vector

    % Create the 3D quiver plot
    quiver3(ax, X, Y, Z, U, V, W, scale, 'Displayname', displayname, 'Color', color, 'LineWidth', 1);
    
    % Set the view angle for better visualization
    view(45, 30);

    % Return axis handles if requested
    if nargout == 2
        varargout{1} = ax.Parent;
        varargout{2} = ax;
    elseif nargout == 1
        varargout{1} = ax;
    end
end
