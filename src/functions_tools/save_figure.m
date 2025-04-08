function [fig, ax] = save_figure(fig, figname, ax)
    % SAVE_FIGURE Saves the current figure with a specified name in EPS format.
    %
    %   [fig, ax] = SAVE_FIGURE(fig, figname, ax)
    %
    % Inputs:
    %   fig     - Handle to the figure to be saved. If not provided, the current figure is used.
    %   figname - Name of the file to save the figure as. If not provided, the current directory name is used.
    %   ax      - Handle to the axes to be saved. If not provided, the current axes are used.
    %
    % Outputs:
    %   fig     - Handle to the figure that was saved.
    %   ax      - Handle to the axes that were saved.
    %
    % Example:
    %   [fig, ax] = save_figure(gcf, 'my_figure', gca);
    %
    % Note:
    %   The function saves the figure in EPS format suitable for publication.
    %   The 'plot_results' directory is created if it does not exist.
    %
    % Author: [Your Name]
    % Date: [Date]
    % Version: 1.0
    % License: [License Information]
    % Contact: [Your Contact Information]

    % Validate input parameters
    if nargin < 1
        fig = gcf; % Use current figure if not provided
    end
    if nargin < 2
        % Use the current directory name as the figure name
        dirname = pwd;
        [~, figname, ~] = fileparts(dirname);
        figname = strcat(figname, '.eps');
    end
    if nargin < 3
        ax = gca; % Use current axes if not provided
    end

    % Create 'plot_results' directory if it doesn't exist
    outputDir = 'plot_results';
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    % Prepare the figure name by replacing special characters
    figname = strrep(figname, ' ', '');
    figname = strrep(figname, '*', '_multi_');
    figname = strrep(figname, '-', '_m_');
    figname = strrep(figname, '+', '_p_');
    figname = strrep(figname, '/', '_divid_');
    figname = strrep(figname, '~', '_power_');

    % Set figure and axes properties for publication quality
    set(ax, 'Color', 'none');
    set(fig, 'Color', 'none');
    set(fig, 'InvertHardCopy', 'off');

    % Save the figure as an EPS file in the 'plot_results' directory
    outputFile = fullfile(outputDir, figname);
    print(fig, outputFile, '-depsc', '-tiff', '-cmyk');

    % Restore figure and axes properties
    set(ax, 'Color', 'white');
    set(fig, 'Color', 'white');
end
