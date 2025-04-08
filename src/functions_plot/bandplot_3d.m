function figure_k = bandplot_3d(bandlist, Ecut, kplane, contourf_mode, wt_mode)
    % bandplot_3d - Visualizes a 3D band structure.
    % 
    % Usage: figure_k = bandplot_3d(bandlist, Ecut, kplane, contourf_mode, wt_mode)
    % 
    % Inputs:
    % bandlist      - List of band indices to plot
    % Ecut          - Energy cutoff for the plot (min, max)
    % kplane        - k-point plane for plotting (e.g., [1, 2])
    % contourf_mode - Whether to use filled contour plots (1 for yes)
    % wt_mode       - Mode of operation (0 = standard, 1 = weighted, 2 = custom data)
    
    % Default argument values
    arguments
        bandlist = [1, 2];
        Ecut = [-1.5, 1.5];
        kplane = [1, 2];
        contourf_mode = 0;
        wt_mode = 0;
    end

    % Create figure for the plot
    figure_k = figure();
    
    % Create axes for the plot
    axes1 = axes('Parent', figure_k);
    hold(axes1, 'on');
    set(axes1, 'FontName', 'Helvetica', 'FontSize', 30);

    % Standard mode (wt_mode == 0)
    if wt_mode == 0
        % Load mesh data
        kx_mesh = load('KX.grd');
        ky_mesh = load('KY.grd');
        
        % Initialize mesh data array
        BM_mesh = cell(1, length(bandlist));
        
        % Load band data for each band in the list
        for i = 1:length(bandlist)
            temp_name = strcat('BAND.B', num2str(bandlist(i)), '.grd');
            BM_mesh{i} = load(temp_name);
        end
        
        % Plot the bands as surface plots
        for i = 1:length(bandlist)
            surf(kx_mesh, ky_mesh, BM_mesh{i}, 'Parent', axes1, 'FaceColor', 'interp', ...
                'EdgeColor', 'none', 'FaceAlpha', 0.6);
            hold on;
        end
        
        % If filled contour mode is enabled, plot the contours
        if contourf_mode == 1
            for i = 1:length(bandlist)
                contourf(kx_mesh, ky_mesh, BM_mesh{i}, 4);
                hold on;
            end
        end
        
        % Axis labels
        xlabel('$\it{k}_{x} (\textrm{1/\AA})$', 'Interpreter', 'latex', 'Fontname', 'Times New Roman');
        ylabel('$\it{k}_{y} (\textrm{1/\AA})$', 'Interpreter', 'latex', 'Fontname', 'Times New Roman');
        zlabel('Energy (eV)');
        axis image;
        axis vis3d;
        shading interp;
        colormap(hsv);
        set(gca, 'FontSize', 20);
        zlim(Ecut);
    
    % Weighted data mode (wt_mode == 1)
    elseif wt_mode == 1
        % Load weighted plane data
        bulkek_plane_matlab = load('bulkek_plane-matlab.dat');
        [~, width] = size(bulkek_plane_matlab);
        
        % Extract k-point and band data
        kx_list = bulkek_plane_matlab(:, 1);
        ky_list = bulkek_plane_matlab(:, 2);
        kz_list = bulkek_plane_matlab(:, 3);
        
        % Number of k-points
        num_k = sqrt(length(kx_list));
        
        % Reshape bulkek data into matrices for each band
        BM_mesh{1} = reshape(bulkek_plane_matlab(:, 7), num_k, num_k);
        BM_mesh{2} = reshape(bulkek_plane_matlab(:, 8), num_k, num_k);
        if width > 8
            BM_mesh{3} = reshape(bulkek_plane_matlab(:, 9), num_k, num_k);
            BM_mesh{4} = reshape(bulkek_plane_matlab(:, 10), num_k, num_k);
        end
        
        % Reshape k-point data
        kx_mesh = reshape(bulkek_plane_matlab(:, kplane(1)), num_k, num_k);
        ky_mesh = reshape(bulkek_plane_matlab(:, kplane(2)), num_k, num_k);
        
        % Plot the bands as surface plots
        for i = bandlist
            surf(kx_mesh, ky_mesh, BM_mesh{i}, 'Parent', axes1, 'FaceColor', 'interp', ...
                'EdgeColor', 'none', 'FaceAlpha', 0.6);
            hold on;
        end
        
        % If filled contour mode is enabled, plot the contours
        if contourf_mode == 1
            for i = bandlist
                contourf(kx_mesh, ky_mesh, BM_mesh{i}, 4);
                hold on;
            end
        end
        
        % Axis labels
        xlabel('k_x ', 'Fontname', 'Helvetica');
        ylabel('k_y', 'Fontname', 'Helvetica');
        zlabel('Energy (eV)');
        axis image;
        axis vis3d;
        shading interp;
        colormap(hsv);
        set(gca, 'FontSize', 30);
        zlim(Ecut);
        view(axes1, [-72, 17]);
    
    % Custom data mode (wt_mode == 2)
    elseif wt_mode == 2
        disp('For MATLAB 3D plot...');
        
        % Reshape the custom data
        num_k = sqrt(length(kplane));
        [nband, ~] = size(bandlist);
        kx_mesh = reshape(kplane(:, 1), num_k, num_k);
        ky_mesh = reshape(kplane(:, 2), num_k, num_k);
        
        % Plot each band as a surface
        for Ei = 1:nband
            E_mesh = reshape(bandlist(Ei, :), num_k, num_k);
            surf(kx_mesh, ky_mesh, E_mesh, 'Parent', axes1, 'FaceColor', 'interp', ...
                'EdgeColor', 'none', 'FaceAlpha', 0.6, 'DisplayName', num2str(Ei));
            hold on;
        end
        
        % If filled contour mode is enabled, plot the contours
        if contourf_mode == 1
            for Ei = 1:nband
                E_mesh = reshape(bandlist(Ei, :), num_k, num_k);
                contourf(kx_mesh, ky_mesh, E_mesh, 4);
                hold on;
            end
        end
        
        % Axis labels
        xlabel('k_x ', 'Fontname', 'Helvetica');
        ylabel('k_y', 'Fontname', 'Helvetica');
        zlabel('Energy (eV)');
        axis image;
        axis vis3d;
        shading interp;
        colormap(hsv);
        set(gca, 'FontSize', 30);
        zlim(Ecut);
        view(axes1, [-72, 17]);
    end
end
