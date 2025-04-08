function ax = ShowSurfIntegral(Gk, klist_r, dk1L, dk2L, dSumL, options)
    % SHOWSURFINTEGRAL Visualizes the surface integral for a given set of inputs.
    %   This function generates three subplots:
    %   1. A plot of the integral patch in the Brillouin zone.
    %   2. A plot of the `dBerryphase` vs. k-path.
    %   3. A cumulative plot of the Berry phase sum and Chern number.
    %
    %   Inputs:
    %       Gk      : The reciprocal lattice vectors (if available).
    %       klist_r : List of k-points.
    %       dk1L    : Vector of first differential values.
    %       dk2L    : Vector of second differential values.
    %       dSumL   : List of differential sums.
    %       options : Struct containing the following optional fields:
    %           - dName       : Name for the differential plot (default: 'dBerryphase').
    %           - SumName     : Name for the cumulative sum plot (default: 'Berryphase').
    %           - FinalName   : Prefix for the final title (default: 'ChernNumber = ').
    %           - oneshot     : If true, prevents live updates (default: true).
    %           - view        : View angles for the first subplot (default: [0, 90]).
    
    arguments
        Gk = [];
        klist_r = [];
        dk1L = [];
        dk2L = [];
        dSumL = [];
        options.dName = 'dBerryphase';
        options.SumName = 'Berryphase';
        options.FinalName = 'ChernNumber = ';
        options.oneshot = true;
        options.view = [0, 90];
    end

    % Create figure and axes for subplots
    fig = figure('PaperType', 'a4letter', 'PaperSize', [16 8], 'Color', 'white', ...
                 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.6]);

    % Subplot 1: Integral Patch in Brillouin Zone
    ax0 = subplot(1, 3, 1, 'LineWidth', 1.5, 'FontSize', 24, 'FontName', "Helvetica", 'Parent', fig);
    hold(ax0, 'all');
    [fig, ax0] = TBkit.BZplot(Gk, 'Gk', true, 'ax', ax0, 'fig', fig, 'color', 'r', 'alpha', 0.3);
    title(ax0, 'integral patch');
    axis(ax0, 'equal');
    view(ax0, options.view(1), options.view(2));

    % Subplot 2: Berry Phase (dName)
    ax1 = subplot(1, 3, 2, 'LineWidth', 1.5, 'FontSize', 24, 'FontName', "Helvetica", 'Parent', fig);
    hold(ax1, 'all');
    xlabel(ax1, 'kpath num');
    ylabel(ax1, options.dName);

    % Subplot 3: Cumulative Berry Phase (SumName)
    ax2 = subplot(1, 3, 3, 'LineWidth', 1.5, 'FontSize', 24, 'FontName', "Helvetica", 'Parent', fig);
    hold(ax2, 'all');
    xlabel(ax2, 'kpath num');
    ylabel(ax2, options.SumName);

    % Prepare the K-points and integrate
    kn = size(klist_r, 1);
    dk1dk2L = dk1L + dk2L;  % Sum of dk1L and dk2L for each k-point
    KsmallL = zeros(kn, 3, 4);
    KsmallL(:, :, 1) = klist_r;
    KsmallL(:, :, 2) = klist_r + dk1L;
    KsmallL(:, :, 3) = klist_r + dk2L;
    KsmallL(:, :, 4) = klist_r + dk1dk2L;

    SUM = 0;
    dSUM_old = 0;
    dSumL = real(dSumL);  % Ensure real values for dSumL

    % Loop over k-points to calculate and update plots
    for ki = 1:kn
        K = klist_r(ki, :);
        dSUM = dSumL(ki);
        
        % Plot the integral patch in the first subplot
        Ksmall = [KsmallL(ki, :, 1); KsmallL(ki, :, 2); KsmallL(ki, :, 4); KsmallL(ki, :, 3)];
        patch(ax0, Ksmall(:, 1), Ksmall(:, 2), Ksmall(:, 3), 'g', 'EdgeColor', 'none');

        % Plot the differential Berry phase (dBerryphase) in the second subplot
        xlabel(ax1, sprintf('%d,%d,%d', K(1), K(2), K(3)));
        area(ax1, [ki - 1, ki], [dSUM_old, dSUM]);
        dSUM_old = dSUM;

        % Update the cumulative Berry phase (Berryphase) in the third subplot
        SUM = SUM + dSUM;
        stem(ax2, ki, SUM);

        % Optionally update the plot live
        if ~options.oneshot
            drawnow;
        end
    end

    % Final title with the calculated Chern number
    title(ax2, [options.FinalName, num2str(SUM)]);
    axis(ax0, 'equal');

    % Return the axes handles
    ax = [ax0, ax1, ax2];
end
