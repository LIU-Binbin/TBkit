function varargout = dos_plot_set(EIGENCAR_list, WEIGHTCAR_list, Name_list, Ecut, titlestring, options)
    % DOS_PLOT_SET Generates a Density of States (DOS) plot from the provided data.
    %   Inputs:
    %       EIGENCAR_list: Matrix containing energy levels (eigenvalues).
    %       WEIGHTCAR_list: Matrix containing corresponding DOS weights.
    %       Name_list: List of names (labels) for the datasets.
    %       Ecut: Energy range for the plot (in eV) [y-axis limits].
    %       titlestring: Title for the plot.
    %       options: Struct containing options, including:
    %           - cmap: Colormap to use for the plot (default: jet).
    %
    %   Outputs:
    %       If nargout >= 1, returns axes (and optionally figure).
    % 

    % Default colormap setup (if not provided in options)
    arguments
        EIGENCAR_list
        WEIGHTCAR_list
        Name_list
        Ecut
        titlestring
        options.cmap = @jet;
    end

    % Initialize the figure
    [~, ax] = Figs('Position', [0.2, 0.2, 0.2, 0.6]);
    
    % Number of datasets
    Nname = size(Name_list, 1);

    % Select colormap based on the number of datasets
    cmap = options.cmap(Nname);

    % Plot each dataset with the corresponding color
    for i = 1:Nname
        plot(ax, WEIGHTCAR_list(:, i), EIGENCAR_list(:, i), ...
            'Color', cmap(i, :), 'LineWidth', 2.0, 'DisplayName', Name_list(i,:));
    end

    % Set plot properties
    ax.YLim = Ecut;  % Y-axis limits based on provided energy cut
    ylabel('E - E_f (eV)');  % Label for y-axis
    legend();  % Display legend
    title(titlestring);  % Title for the plot
    
    % Return axes or figure if requested
    if nargout == 2
        varargout{1} = ax.Parent;
        varargout{2} = ax;
    elseif nargout == 1
        varargout{1} = ax;
    end
end
