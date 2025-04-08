function varargout = create_figure(M, N, propArgs, FigureOption)
    % CREATE_FIGURE Creates a figure with customizable properties and size.
    %   This function generates a figure and optionally sets axes and figure properties.
    %   The function allows specifying the number of rows (M) and columns (N) for subplots
    %   as well as various customizable properties such as figure position and others.
    %
    % Inputs:
    %   M              - Number of rows for subplots (default: 1).
    %   N              - Number of columns for subplots (default: 1).
    %   propArgs       - A structure containing additional properties for figure/axes creation.
    %                    https://ww2.mathworks.cn/help/matlab/matlab_prog/function-argument-validation-1.html#mw_1b62b6d6-a445-4c55-a9b9-9c70becfdbe6      
    %       test       - Placeholder for any additional options in propArgs (optional).
    %   FigureOption   - Options for figure customization.
    %       Position   - Position of the figure on the screen [default: empty, meaning no position change].
    %
    % Outputs:
    %   varargout      - Depending on the number of output arguments requested, 
    %                    returns the figure handle, axes handles, or the full Fig structure.

    % Argument validation and default value assignment.
    arguments
        M = 1;                     % Default to 1 row.
        N = 1;                     % Default to 1 column.
        propArgs.test = [];        % Placeholder for additional arguments (can be extended).
        FigureOption.Position = []; % Default empty position (no change).
    end
    
    % Convert property arguments structure to cell array for flexible usage.
    propArgsCell = namedargs2cell(propArgs);  
    
    % Create a new figure with subplots defined by M (rows) and N (columns).
    Fig = Figs(M, N);  % Figs is assumed to be a function that creates a figure and axes.

    %-------- Initialization Section --------
    % If a position for the figure is specified, apply it.
    if ~isempty(FigureOption.Position)
        set(Fig.fig, 'Position', FigureOption.Position);  % Set figure position if defined.
    end
    
    %-------- Return Section --------
    % Return the figure and axes depending on the number of output arguments requested.
    if nargout == 2
        varargout{1} = Fig.fig;     % Return the figure handle.
        varargout{2} = Fig.axes;    % Return the axes handles.
    elseif nargout == 1
        varargout{1} = Fig;         % Return the full Fig object (figure and axes).
    end
end
    
