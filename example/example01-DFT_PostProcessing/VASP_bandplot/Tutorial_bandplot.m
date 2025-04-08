%% MATLAB Band Structure Plotting Tutorial
% This tutorial provides step-by-step instructions for plotting the band structure of materials
% using the `bandplot` function in MATLAB. Ensure that the necessary input files are in the current directory
% before proceeding.

%% Step 1: Load Required Files
% Ensure that the following files are available in your current working directory:
% - POSCAR: Contains atomic structure information
% - KPOINTS: Defines the k-point sampling for the band structure
% - EIGENVAL: Contains eigenvalue data for the bands
% - DOSCAR: (Optional) Contains density of states data

% Use the 'ls' command to list the files in the current directory:
ls 

%% Step 2: Default Band Plot
% To create a basic band structure plot using default parameters, you can simply call:
%   bandplot
% or 
%   bandplot()
% This will generate a plot using default settings.

bandplot; % Calls the bandplot function using default parameters.

%% Step 3: User-Defined Band Plot
% Before creating a custom band plot, load the eigenvalue data using the `EIGENVAL_read()` function.
% This will produce a structured matrix known as `EIGENCAR`, which contains the eigenvalue data.

EIGENCAR = EIGENVAL_read(); % Load eigenvalue data from the output of a VASP calculation.

%% Step 4: Define Parameters for the Band Plot
% The full syntax for calling the `bandplot` function with various parameters is as follows:
%   [fig, ax] = bandplot(EIGENCAR, Ecut, titlestring, color, klist_l, kpoints_l, kpoints_name, fontname, fig, ax)

% For this example:
% - Set the energy cut range: Ecut = [-6, 6]
% - Define the title for the plot as 'Band Structure CuI w SOC'
% - Specify color for the bands as 'blue' ('b')

% Call the `bandplot` function with the defined parameters:
bandplot(EIGENCAR, [-6, 6], 'title', 'Band Structure CuI w SOC', 'Color', 'b');

%% Step 5: Save Plot Data
% To obtain the figure handle and the plot data created by `bandplot`, you can store the outputs as follows:
[fig, ax, data_plot] = bandplot(); % Get figure, axis, and plot data for further use.

% After creating the plot, you can export it to different formats. For example, to save the current figure as a PNG:
export_fig(gcf, 'plot_results/CuI_band.png'); % Save figure as PNG format (you can also save as pdf, eps, etc.)

% To save the figure in SVG format with high resolution, use the following command:
set(gcf, 'Renderer', 'Painters'); % Set the renderer to 'Painters' for better SVG output.
print(gcf, '-dsvg', '-r600', 'plot_results/CuI_band.svg'); % Save the figure as SVG with 600 DPI resolution.

%% Step 6: Multi-Bandplot  Plotting
% In this step, we will create a multi-band structures plot to visualize multiple sets of eigenvalue data on the same axes.
% This is useful for comparing the band structures of different materials or the same material under different conditions.

% First, we will assign the original eigenvalue data to a variable.
EIGENCAR1 = EIGENCAR; % Set the first set of eigenvalues.

% Next, we will create a second set of eigenvalues by simply shifting the original data upwards by 0.03 eV.
EIGENCAR2 = EIGENCAR + 0.03; % Set the second set of eigenvalues, slightly altered for the comparison.

% Now we will plot the first band structure using the `bandplot` function.
ax = bandplot(EIGENCAR1, [-6, 6], ...
    'title', 'Band Structures', 'Color', 'b', ...  % Set the title and color to blue
    'LineSpec', '-'); % Use solid lines for the first plot

% Next, we will plot the second band structure on the same axes created by the first call to `bandplot`.
% We reuse the axis handle (`ax`) to ensure that both plots appear on the same figure.
ax = bandplot(EIGENCAR2, [-6, 6], 'ax', ax, ... % Use the axes handle from the first plot
    'title', 'Band Structures', 'Color', 'r', ... % Set the title and color to red for contrast
    'LineSpec', '--', ... % Use dashed lines for the second plot
    'LineWidth', 0.5); % Adjust the line width for visibility

% At this point, both band structures will be displayed on the same figure, allowing for easy comparison.
% The blue solid lines represent the first set of eigenvalues, while the red dashed lines represent the shifted eigenvalues.

