function varargout = waveplot(orb_list, WaveFunc, options)
    % WAVEPLOT Plots the wavefunction with orbital positions and magnitudes.
    %   The function displays the orbital positions along with their wavefunction
    %   magnitudes, allowing customization of colors, sizes, and the option for the 
    %   POSCAR file input for the crystal structure.
    %
    %   Inputs:
    %       orb_list : Nx3 matrix of orbital positions.
    %       WaveFunc : NxM matrix of wavefunctions (M wavefunctions, N orbitals).
    %       options  : A structure with the following optional fields:
    %           - ax           : Axes handle for the plot (optional).
    %           - Rm           : Unit cell (optional).
    %           - POSCAR       : POSCAR filename (optional, default 'POSCAR').
    %           - WaveMin      : Minimum wavefunction value for plotting (default 1e-3).
    %           - WaveColor    : Color of the wavefunction scatter plot (default 'r').
    %           - WaveSize     : Size of wavefunction markers (default 1).
    %           - OrbColor     : Color of the orbital scatter plot (default 'k').
    %           - OrbSize      : Size of orbital markers (default 1).
    %
    %   Outputs:
    %       varargout  : Depending on the number of output arguments, returns axes handles.

    arguments
        orb_list = [];                      % Orbital list (Nx3 matrix)
        WaveFunc = zeros(size(orb_list,1), 1); % Wavefunction (NxM matrix)
        options.ax = handle([]);            % Axes handle (optional)
        options.Rm = [];                    % Unit cell (optional)
        options.POSCAR = 'POSCAR';          % POSCAR file (optional)
        options.WaveMin = 1e-3;             % Minimum wavefunction magnitude to plot (default 1e-3)
        options.WaveColor = 'r';            % Color of wavefunction markers (default 'r')
        options.WaveSize = 1;               % Size of wavefunction markers (default 1)
        options.OrbColor = 'k';             % Color of orbital markers (default 'k')
        options.OrbSize = 1;                % Size of orbital markers (default 1)
    end

    % Read POSCAR if Rm is not provided
    if isempty(options.Rm)
        try
            [Rm,~,~,~,~] = POSCAR_read(options.POSCAR); % Try reading POSCAR file
        catch
            filename = input('Please provide a POSCAR file: '); % Prompt user for POSCAR file
            [Rm,~,~,~,~] = POSCAR_readin(filename);  % Read the provided file
        end
    else
        Rm = options.Rm;  % Use the provided unit cell if available
    end

    % Set axes handle
    if isempty(options.ax)
        [~,ax] = Figs(1, 1);  % Create a new figure if no axes handle is provided
    else
        if ishandle(options.ax)
            ax = options.ax;  % Use the provided axes handle if valid
        else
            error('Invalid axes handle provided.');
        end
    end

    % Validate wavefunction and orbital list dimensions
    [Nwave, Nlist] = size(WaveFunc);
    Norb = length(orb_list);
    if Norb ~= Nwave
        error('Orbital list length is not equal to WaveFunc.');
    end

    % Prepare data for plotting
    WFplot_list = zeros(Norb, 4);
    WFplot_list(:, 1:3) = orb_list * Rm;  % Orbital positions scaled by Rm
    WFplot_list(:, 4) = sum(WaveFunc .* conj(WaveFunc), 2);  % Wavefunction magnitude squared

    % Create mask for displaying wavefunction values above threshold
    label_list = WFplot_list(:, 4) > options.WaveMin;  % Threshold filtering
    label_list2 = ~label_list;  % For orbital markers

    % Plot wavefunction values
    scatter3(ax, WFplot_list(label_list, 1), WFplot_list(label_list, 2), WFplot_list(label_list, 3), ...
        WFplot_list(label_list, 4) * 1000 * options.WaveSize, 'filled', 'MarkerFaceColor', options.WaveColor);

    hold(ax, 'on');
    grid(ax, 'off');

    % Plot orbital positions with minimal wavefunction
    plot3(ax, WFplot_list(label_list2, 1), WFplot_list(label_list2, 2), WFplot_list(label_list2, 3), ...
        'ko', 'MarkerSize', options.OrbSize, 'MarkerFaceColor', options.OrbColor);

    % Set view and axis properties
    view(0, 90);  % Top-down view
    axis(ax, 'equal');
    axis(ax, 'off');

    % Return axes handles if requested
    if nargout == 2
        varargout{1} = ax.Parent;
        varargout{2} = ax;
    elseif nargout == 1
        varargout{1} = ax;
    end
end
