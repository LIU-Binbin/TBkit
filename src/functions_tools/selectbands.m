function [num_bands, Nmin, Nmax,win_info] = selectbands(EIGENCAR, Rm, mode,options)
    % selectbands - Select and plot energy bands or check energy window
    %
    % Syntax:
    %   [num_bands, Nmin, Nmax] = selectbands(EIGENCAR, Rm, mode)
    %
    % Input:
    %   EIGENCAR  - The energy bands data (NxM matrix), where N is the number of bands and M is the number of k-points
    %   Rm        - Crystal lattice parameters (3x3 matrix) used for k-path generation
    %   mode      - String specifying the mode of operation ('win_gen' or 'dis_check'). Default is 'win_gen'.
    %
    % Output:
    %   num_bands - Number of selected energy bands (Nmax - Nmin + 1)
    %   Nmin       - The minimum band index selected by the user
    %   Nmax       - The maximum band index selected by the user
    arguments
        EIGENCAR 
        Rm  = eye(3);
        mode {mustBeMember(mode,{'win_gen','dis_check','win_gen&dis_check'})} = 'win_gen'; % Default mode is 'win_gen'
        options.BandIndex = []; 
    end
    
    % Prompt user for Nmin and Nmax input and validate the input
    if isempty(options.BandIndex)
        while true
            Nmin = input('Please input the Nmin (minimum band index):\n');
            Nmax = input('Please input the Nmax (maximum band index):\n');

            % Validate Nmin and Nmax input
            if Nmin < 1 || Nmax > size(EIGENCAR, 1) || Nmin > Nmax
                disp('Invalid input! Ensure that Nmin and Nmax are within valid band indices and Nmin <= Nmax.');
            else
                break; % Exit the loop when valid input is provided
            end
        end
    else
        Nmin = options.BandIndex(1);
        Nmax = options.BandIndex(2);
    end
    
    % Mode 'win_gen' - Generate and plot energy bands
    if strcmp(mode, 'win_gen')
        num_bands = Nmax - Nmin + 1;
        
        % Generate k-path for 3D data
        [klist_r, klist_l, klist_s, kpoints_l, kpoints_name] = kpathgen3D(Rm);
        size_ = [0.2, 0.2, 0.6, 0.6]; % Position of the plot
        fontname = 'Helvetica'; % Font used in the plot
        
        kpathname = kpoints_name;
        klist = klist_l;
        
        % Create the figure for plotting the bands
        figure_k = figure('PaperType', 'a4letter', 'PaperSize', [8 6], 'Color', [1 1 1]);
        axes1 = axes('Parent', figure_k, 'Position', size_, 'LineWidth', 1.5, 'FontSize', 30.0, 'FontName', fontname);
        box(axes1, 'on');
        hold(axes1, 'all');
        
        % Plot all the energy bands
        Nbands = size(EIGENCAR, 1);
        for Ei = 1:Nbands
            plot(klist, EIGENCAR(Ei, :), 'LineWidth', 2.0, 'Color', [0 0 1], 'DisplayName', num2str(Ei));
            hold on;
        end
        
        % Highlight the selected bands
        plot(klist, EIGENCAR(Nmin, :), 'LineWidth', 3.0, 'Color', [1 0 0], 'DisplayName', 'Nmin');
        hold on;
        plot(klist, EIGENCAR(Nmax, :), 'LineWidth', 3.0, 'Color', [1 0 0], 'DisplayName', 'Nmax');
        hold on;
        plot(klist, EIGENCAR(Nmin - 1, :), 'LineWidth', 3.0, 'Color', [1 1 0], 'DisplayName', 'Nmin', 'LineStyle', '--');
        hold on;
        plot(klist, EIGENCAR(Nmax + 1, :), 'LineWidth', 3.0, 'Color', [1 1 0], 'DisplayName', 'Nmax', 'LineStyle', '--');
        hold on;
        
        % Set the axis limits and labels
        xmin = kpoints_l(1);
        xmax = kpoints_l(end);
        ymin = min(min(EIGENCAR(:,:)));
        ymax = max(max(EIGENCAR(:,:)));
        axis([xmin, xmax, ymin, ymax]);
        
        % Replace 'GAMMA' and 'G' with \Gamma for better readability
        kpathname = strrep(kpathname, 'GAMMA', 'G');
        kpathname = strrep(kpathname, 'G', '\Gamma');
        
        % Set plot labels and tick marks
        set(gca, 'XTickLabel', '');
        set(axes1, 'FontName', 'Helvetica', 'FontSize', 30, 'LineWidth', 1.5, 'XTick', kpoints_l, 'XTickLabel', kpathname);
        ylabel('E-E_F (eV)', 'FontWeight', 'bold', 'FontName', fontname);
        
        % Add Fermi level and k-path lines
        line([xmin, xmax], [0, 0], 'LineWidth', 0.5, 'Color', [0 0 0], 'DisplayName', 'Fermi level');
        for i = 1:(length(kpathname) - 2)
            line([kpoints_l(i + 1), kpoints_l(i + 1)], [ymin, ymax], 'LineWidth', 0.2, 'Color', [0 0 0], 'DisplayName', 'K-path');
        end
        
        % Add title
        title('Energy Bands with Selected Range');
        
    % Mode 'dis_check' - Check the energy window and suggest a range
    elseif strcmp(mode, 'dis_check')
        infinite_small = 1e-2; % Small margin for the energy window
        
        % Calculate the energy window for the selected bands
        dis_min = min(EIGENCAR(Nmin, :)) - infinite_small;
        dis_max = max(EIGENCAR(Nmax, :)) + infinite_small;
        
        % Use findbands (assumed to be defined) to get the bands within the window
        [allk_min, anyk_max] = findbands(EIGENCAR, [dis_min, dis_max]);
        
        % Suggested energy window for frozen bands
        foz_min = max(EIGENCAR(Nmin, :)) + infinite_small;
        foz_max = min(EIGENCAR(Nmax, :)) - infinite_small;
        
        % Display recommendations to the user
        disp('Recommended band window:');
        fprintf('Nbands = %d !',anyk_max.nbands);
        fprintf('choose (%d - %d)\n',anyk_max.n1,anyk_max.n2);
        fprintf('dis_win_min = %f\n',dis_min);
        fprintf('dis_win_max = %f\n',dis_max);
        fprintf('dis_froz_min = %f\n',foz_min);
        fprintf('dis_froz_max = %f\n',foz_max);
        fprintf('another choice: \n')
        fprintf('dis_win_min < %f\n',dis_min);
        fprintf('dis_win_max > %f\n',dis_max);
        fprintf('dis_froz_min > %f\n',foz_min);
        fprintf('dis_froz_max < %f\n',foz_max);
        win_info.dis_min = dis_min;
        win_info.dis_max = dis_max;
        win_info.foz_min = foz_min;
        win_info.foz_max = foz_max;
                Nmin = anyk_max.n1;
        Nmax = anyk_max.n2;
        num_bands = anyk_max.nbands;
    elseif strcmp(mode, 'win_gen&dis_check')
        infinite_small = 1e-2; % Small margin for the energy window
        
        % Calculate the energy window for the selected bands
        dis_min = min(EIGENCAR(Nmin, :)) - infinite_small;
        dis_max = max(EIGENCAR(Nmax, :)) + infinite_small;
        
        % Use findbands (assumed to be defined) to get the bands within the window
        [allk_min, anyk_max] = findbands(EIGENCAR, [dis_min, dis_max]);
        
        % Suggested energy window for frozen bands
        foz_min = max(EIGENCAR(Nmin, :)) + infinite_small;
        foz_max = min(EIGENCAR(Nmax, :)) - infinite_small;
        
        % Display recommendations to the user
        disp('Recommended band window:');
        fprintf('Nbands = %d !',anyk_max.nbands);
        fprintf('choose (%d - %d)\n',anyk_max.n1,anyk_max.n2);
        fprintf('dis_win_min = %f\n',dis_min);
        fprintf('dis_win_max = %f\n',dis_max);
        fprintf('dis_froz_min = %f\n',foz_min);
        fprintf('dis_froz_max = %f\n',foz_max);
        fprintf('another choice: \n')
        fprintf('dis_win_min < %f\n',dis_min);
        fprintf('dis_win_max > %f\n',dis_max);
        fprintf('dis_froz_min > %f\n',foz_min);
        fprintf('dis_froz_max < %f\n',foz_max);
        win_info.dis_min = dis_min;
        win_info.dis_max = dis_max;
        win_info.foz_min = foz_min;
        win_info.foz_max = foz_max;
        % Generate k-path for 3D data
        [klist_r, klist_l, klist_s, kpoints_l, kpoints_name] = kpathgen3D(Rm);
        size_ = [0.2, 0.2, 0.6, 0.6]; % Position of the plot
        fontname = 'Helvetica'; % Font used in the plot
        
        kpathname = kpoints_name;
        klist = klist_l;
        
        % Create the figure for plotting the bands
        figure_k = figure('PaperType', 'a4letter', 'PaperSize', [8 6], 'Color', [1 1 1]);
        axes1 = axes('Parent', figure_k, 'Position', size_, 'LineWidth', 1.5, 'FontSize', 30.0, 'FontName', fontname);
        box(axes1, 'on');
        hold(axes1, 'all');
        
        % Plot all the energy bands
        Nbands = size(EIGENCAR, 1);
        for Ei = 1:Nbands
            plot(klist, EIGENCAR(Ei, :), 'LineWidth', 2.0, 'Color', [0 0 1], 'DisplayName', num2str(Ei));
            hold on;
        end
        Nmin = anyk_max.n1;
        Nmax = anyk_max.n2;
        num_bands = anyk_max.nbands;
        % Highlight the selected bands
        plot(klist, EIGENCAR(anyk_max.n1, :), 'LineWidth', 3.0, 'Color', [1 0 0], 'DisplayName', 'DFT:Nmin');
        hold on;
        plot(klist, EIGENCAR(anyk_max.n2, :), 'LineWidth', 3.0, 'Color', [1 0 0], 'DisplayName', 'DFT:Nmax');
        hold on;
        plot(klist, EIGENCAR(anyk_max.n1 - 1, :), 'LineWidth', 3.0, 'Color', [1 1 0], 'DisplayName', 'Boundary', 'LineStyle', '--');
        hold on;
        plot(klist, EIGENCAR(anyk_max.n2 + 1, :), 'LineWidth', 3.0, 'Color', [1 1 0], 'DisplayName', 'Boundary', 'LineStyle', '--');
        hold on;
        
        % Set the axis limits and labels
        xmin = kpoints_l(1);
        xmax = kpoints_l(end);
        ymin = min(min(EIGENCAR(:,:)));
        ymax = max(max(EIGENCAR(:,:)));
        axis([xmin, xmax, ymin, ymax]);
        
        % Replace 'GAMMA' and 'G' with \Gamma for better readability
        kpathname = strrep(kpathname, 'GAMMA', 'G');
        kpathname = strrep(kpathname, 'G', '\Gamma');
        
        % Set plot labels and tick marks
        set(gca, 'XTickLabel', '');
        set(axes1, 'FontName', 'Helvetica', 'FontSize', 30, 'LineWidth', 1.5, 'XTick', kpoints_l, 'XTickLabel', kpathname);
        ylabel('E-E_F (eV)', 'FontWeight', 'bold', 'FontName', fontname);
        
        % Add Fermi level and k-path lines
        line([xmin, xmax], [0, 0], 'LineWidth', 0.5, 'Color', [0 0 0], 'DisplayName', 'Fermi level');
        for i = 1:(length(kpathname) - 2)
            line([kpoints_l(i + 1), kpoints_l(i + 1)], [ymin, ymax], 'LineWidth', 0.2, 'Color', [0 0 0], 'DisplayName', 'K-path');
        end
        
        % Add title
        title('Energy Bands with Selected Range');
    end
end
