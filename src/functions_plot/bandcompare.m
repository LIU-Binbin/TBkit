function ax = bandcompare(EIGENCAR, EIGENCAR_DFT, Ecut, klist_l, kpoints_l, kpoints_name, options)
    % bandcompare - Compare band structures from two different sources and plot them.
    % 
    % This function plots the electronic band structures from two data sets:
    % EIGENCAR and EIGENCAR_DFT. Both sets are compared along the same k-point
    % path, and the corresponding energy values are displayed on the same graph.
    %
    % Arguments:
    % EIGENCAR        - The EIGENCAR data (can be output from a VASP run).
    % EIGENCAR_DFT    - The DFT EIGENCAR data for comparison (default is read from a file).
    % Ecut            - Energy range for the plot (default: [-3, 3]).
    % klist_l         - List of k-points for the plot.
    % kpoints_l       - List of k-point indices.
    % kpoints_name    - Names of the k-points.
    % options         - A structure containing additional plotting options.
    %
    % Returns:
    % ax              - Handle to the axis object containing the plot.

    % Input parsing with default values
    arguments
        EIGENCAR;
        EIGENCAR_DFT = EIGENVAL_read();  % Default function to read EIGENVAL file
        Ecut double = [-3, 3];           % Default energy cutoff range for plotting
        klist_l double = [];             % Optional k-point list (empty by default)
        kpoints_l double = [];           % Optional k-point indices (empty by default)
        kpoints_name string = [];        % Optional k-point names (empty by default)
        options.ax = handle([]);         % Axis handle, empty by default
        options.fontname = 'Helvetica';  % Font for labels
        options.KPOINTS = 'KPOINTS';     % KPOINTS file (if relevant)
        options.POSCAR = 'POSCAR';       % POSCAR file (if relevant)
        options.Color = @jet;            % Color map for plotting (default: jet)
        options.title = '';              % Title for the plot (default: empty)
        options.klist_l = [];            % klist_l (overrides if provided)
        options.kpoints_l = [];          % kpoints_l (overrides if provided)
        options.kpoints_name = [];       % kpoints_name (overrides if provided)
        options.xlabel = '';             % Label for x-axis
        options.ylabel = 'E (eV)';       % Label for y-axis
        options.LineSpec = '-';          % Line style (default: solid line)
        options.LineWidth = 1;           % Line width (default: 1)
        options.MarkerSize = 3;          % Marker size (default: 3)
        options.MarkerEdgeColor = [];    % Marker edge color (default: [])
        options.MarkerFaceColor = [];    % Marker face color (default: [])
    end

    % Default plotting settings (ensure they're set in case options are empty)
    options.LineSpec = '-';
    options.Color = 'b';  % Default color for EIGENCAR_DFT plot

    % Convert options to a cell array for property name-value pairs
    propertyCell1 = namedargs2cell(options);

    % Plot the DFT data first
    options.ax = bandplot(EIGENCAR_DFT, Ecut, klist_l, kpoints_l, kpoints_name, propertyCell1{:});
    
    % Change settings for EIGENCAR data plot (red color and dot style)
    options.LineSpec = '.';  % Change to dotted line style for the second plot
    options.Color = 'r';     % Change to red color for EIGENCAR plot

    % Convert options to a cell array for the second plot
    propertyCell2 = namedargs2cell(options);

    % Plot the EIGENCAR data on top of the DFT data
    ax = bandplot(EIGENCAR, Ecut, klist_l, kpoints_l, kpoints_name, propertyCell2{:});
    
    % Return the axis handle
end
