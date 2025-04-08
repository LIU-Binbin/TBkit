function uiwaveplot(orb_list, WAVECAR_disk, EIGENCAR_disk, norb, Rm, options)
    % UIWAVEPLOT Interactively explore and visualize wavefunction data.
    % This function creates a GUI that allows users to modify parameters for 
    % plotting wavefunctions and explore the energy bands interactively.
    %
    % Inputs:
    %   orb_list        - List of orbital coordinates for plotting.
    %   WAVECAR_disk    - Data from the WAVECAR file (wavefunction data).
    %   EIGENCAR_disk   - Data from the EIGENCAR file (energy bands).
    %   norb            - Number of orbitals (band structure size).
    %   Rm              - Lattice matrix (from POSCAR or similar).
    %   options         - A struct containing options such as:
    %                     Ecut: Energy cut for plot limits.
    %                     Scale: Scaling factor for wavefunction visualization.
    %                     NumBands: Number of bands to display.
    %                     WaveMin: Minimum value for wavefunction visualization.
    %                     selectbands: Initial band selection.
    %                     Position: Window position and size.
    
    arguments
        orb_list
        WAVECAR_disk
        EIGENCAR_disk
        norb
        Rm
        options.Ecut = [-6, 1];        % Energy cut for plotting.
        options.Scale = 1;             % Scaling factor for wavefunctions.
        options.NumBands = 1;         % Number of bands to display.
        options.WaveMin = 1e-3;       % Minimum wavefunction magnitude.
        options.selectbands = [];     % Initial selection of bands.
        options.Position = [0, 100, 1420, 980];  % Window position and size.
    end
    
    % Create the figure window for the interactive plot.
    fig = uifigure('Position', options.Position);
    global selectband scatter_handle;
    
    % Set initial band selection if not provided.
    if isempty(options.selectbands)
        selectband = ceil(norb / 2);
    else
        selectband = options.selectbands;
    end
    
    % Extract parameter values from options struct.
    Scale = options.Scale;
    NumBands = options.NumBands;
    WaveMin = options.WaveMin;
    
    % Set up the layout for the GUI using uigridlayout.
    gl = uigridlayout(fig, [4, 10]);
    gl.RowHeight = {22, 22, 22, '1x'};
    gl.ColumnWidth = {105, 105, 105, 105, '1x'};
    
    % Position axes for plotting energy band data and wavefunction.
    ax1 = uiaxes(gl);
    ax1.Layout.Row = 4;
    ax1.Layout.Column = [1 4];
    
    ax2 = uiaxes(gl);
    ax2.Layout.Row = [1 4];
    ax2.Layout.Column = [5 10];
    
    % Energy range for the plot.
    Ecut = options.Ecut;
    
    % Create UI components (buttons and editable fields) for interaction.
    
    % Button to navigate to the previous energy band.
    btn1 = uibutton(gl, "ButtonPushedFcn", @(src, event) plotButtonPushed1(ax1, ax2, Ecut));
    btn1.Text = 'Previous';
    btn1.Layout.Row = 1;
    btn1.Layout.Column = [1 2];
    
    % Button to navigate to the next energy band.
    btn2 = uibutton(gl, "ButtonPushedFcn", @(src, event) plotButtonPushed2(ax1, ax2, Ecut));
    btn2.Text = 'Next';
    btn2.Layout.Row = 1;
    btn2.Layout.Column = [3 4];
    
    % Editable field for wavefunction scale.
    Field1 = uieditfield(gl, "numeric", ...
        "Value", Scale, ...
        "ValueDisplayFormat", "Scale: %.2f", ...
        "Limits", [0, inf], ...
        "ValueChangedFcn", @(src, event) editFieldValueChanged1(src, event));
    Field1.Layout.Row = 2;
    Field1.Layout.Column = [1 2];
    
    % Editable field for number of bands to display.
    Field2 = uieditfield(gl, "numeric", ...
        "Value", NumBands, ...
        "ValueDisplayFormat", "NumBands: %d", ...
        "Limits", [0, inf], ...
        "ValueChangedFcn", @(src, event) editFieldValueChanged2(src, event));
    Field2.Layout.Row = 2;
    Field2.Layout.Column = [3 4];
    
    % Editable field for minimum wavefunction value.
    Field3 = uieditfield(gl, "numeric", ...
        "Value", WaveMin, ...
        "ValueDisplayFormat", "WaveMin: %.3e", ...
        "Limits", [0, inf], ...
        "ValueChangedFcn", @(src, event) editFieldValueChanged3(src, event));
    Field3.Layout.Row = 3;
    Field3.Layout.Column = [1 4];
    
    % Initial plot of the energy bands.
    plot(ax1, 1:norb, EIGENCAR_disk, '-o');
    hold(ax1, 'on');
    title(ax1, ['Select State: ', num2str(selectband), ' Energy: ', num2str(EIGENCAR_disk(selectband)), ' eV']);
    ylim(ax1, Ecut);
    
    % Plot the wavefunction for the selected band.
    plotAx(ax1, ax2, orb_list, WAVECAR_disk, EIGENCAR_disk, norb, selectband, Rm, options.Ecut);
    
    % Close request function to confirm before closing the window.
    fig.CloseRequestFcn = @(src, event) my_closereq(src);
    fig.Name = 'Select Wave Plot';
    
    % Callback functions for the buttons and input fields.
    
    % When the "Previous" button is pressed, navigate to the previous band.
    function plotButtonPushed1(ax1, ax2, Ecut)
        selectband = selectband - 1;
        plotAx(ax1, ax2, orb_list, WAVECAR_disk, EIGENCAR_disk, norb, selectband, Rm, Ecut);
    end
    
    % When the "Next" button is pressed, navigate to the next band.
    function plotButtonPushed2(ax1, ax2, Ecut)
        selectband = selectband + 1;
        plotAx(ax1, ax2, orb_list, WAVECAR_disk, EIGENCAR_disk, norb, selectband, Rm, Ecut);
    end
    
    % When the Scale field value is changed, update the plot.
    function editFieldValueChanged1(src, event)
        Scale = src.Value;
        plotAx(ax1, ax2, orb_list, WAVECAR_disk, EIGENCAR_disk, norb, selectband, Rm, Ecut);
    end
    
    % When the NumBands field value is changed, update the plot.
    function editFieldValueChanged2(src, event)
        NumBands = src.Value;
        plotAx(ax1, ax2, orb_list, WAVECAR_disk, EIGENCAR_disk, norb, selectband, Rm, Ecut);
    end
    
    % When the WaveMin field value is changed, update the plot.
    function editFieldValueChanged3(src, event)
        WaveMin = src.Value;
        plotAx(ax1, ax2, orb_list, WAVECAR_disk, EIGENCAR_disk, norb, selectband, Rm, Ecut);
    end
    
    % Function to update the plots based on selected band and parameters.
    function plotAx(ax1, ax2, orb_list, WAVECAR_disk, EIGENCAR_disk, norb, selectband, Rm, Ecut)
        % Adjust selected band range based on NumBands.
        if NumBands == 1
            selectband = selectband;
        elseif mod(NumBands, 2) == 0
            selectband = (selectband - floor(NumBands / 2) + 1):(selectband + floor(NumBands / 2));
        else
            selectband = selectband - floor(NumBands / 2):(selectband + floor(NumBands / 2));
        end
        
        % Extract wavefunction data for the selected bands.
        WaveFunc = WAVECAR_disk(:, selectband);
        
        % Clear previous scatter plot if it exists.
        try
            delete(scatter_handle);
        catch
        end
        
        % Plot energy bands on ax1.
        scatter_handle = scatter(ax1, selectband, EIGENCAR_disk(selectband), 20, 'red', 'filled', 'h');
        
        % Plot wavefunction data on ax2.
        waveplot(orb_list, WaveFunc, 'Rm', Rm, 'ax', ax2, 'OrbSize', 3, 'WaveSize', Scale, 'WaveMin', WaveMin);
        view(ax2, 0, 90);
        axis(ax2, 'equal');
        hold(ax2, 'off');
    end
    
    % Function to handle the figure close event with confirmation.
    function my_closereq(fig)
        selection = uiconfirm(fig, 'Close the figure window?', 'Confirmation');
        switch selection
            case 'OK'
                delete(fig);
            case 'Cancel'
                return;
        end
    end
end
