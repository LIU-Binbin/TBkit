function gnu_gen(title, E_cut, POSCAR_file, KPOINTS_file, fermi, mode)
% gnu_gen - Generates a GNUplot script for band structure plotting.
%
% Usage:
%   gnu_gen(title, E_cut, POSCAR_file, KPOINTS_file, fermi, mode)
%
% Inputs:
%   title        - Title for the output plot and script file (default: 'BAND').
%   E_cut        - Energy range for the plot (default: [-3, 3]).
%   POSCAR_file  - Name of the POSCAR file (default: 'POSCAR').
%   KPOINTS_file - Name of the KPOINTS file (default: 'KPOINTS').
%   fermi        - Fermi energy (either as a numeric value or a filename, default: 'DOSCAR').
%   mode         - Mode for the calculations (e.g., 'vasp', 'qe', default: 'vasp').
%
% This function generates a `.gnu` script file for plotting a band structure using GNUplot.

    %% Handle input arguments
    % Set default values for missing arguments
    if nargin < 1, title = 'BAND'; end
    if nargin < 2, E_cut = [-3, 3]; end
    if nargin < 3, POSCAR_file = 'POSCAR'; end
    if nargin < 4, KPOINTS_file = 'KPOINTS'; end
    if nargin < 5, fermi = 'DOSCAR'; end
    if nargin < 6, mode = 'vasp'; end
    
    %% Process Fermi energy
    if strcmp(class(fermi), 'double')
        % Fermi is already a numeric value
    elseif ischar(fermi) || isa(fermi, 'string')
        % Fermi is given as a filename or string; get the value from the file
        fermi = GetFermi(mode, fermi);
    else
        % Invalid Fermi input, set it to 0
        fermi = 0;
    end
    
    %% Read POSCAR and KPOINTS files
    % Read the POSCAR file for lattice information
    [Rm, ~, ~, ~] = POSCAR_readin(POSCAR_file, mode);
    
    % Read the KPOINTS file for k-point information
    [kpoints, nodes, kpoints_name] = KPOINTS_read(KPOINTS_file, mode);
    
    % Generate the k-path in 3D space based on POSCAR and KPOINTS data
    [~, ~, ~, kpoints_l, kpoints_name] = kpathgen3D(Rm, kpoints, nodes, kpoints_name);

    %% Generate the GNUplot script
    % Open a new file for the GNUplot script
    fid = fopen([title, '.gnu'], 'w');
    
    % Write initial setup commands for the plot
    fprintf(fid, 'set terminal postscript enhanced color font ",20"\n');
    fprintf(fid, 'set palette defined (-1 "blue", 0 "grey", 1 "red")\n');
    fprintf(fid, 'set colorbox vertical user origin .1,.35 size .02,.4\n');
    fprintf(fid, 'set output "%s.eps"\n', title);
    fprintf(fid, 'set style data linespoints\n');
    fprintf(fid, 'unset ztics\n');
    fprintf(fid, 'unset key\n');
    fprintf(fid, 'set pointsize 0.8\n');
    fprintf(fid, 'set view 0,0\n');

    % Set axis properties
    fprintf(fid, 'set xtics font ",24"\n');
    fprintf(fid, 'set ytics font ",24"\n');
    fprintf(fid, 'set ylabel font ",24"\n');
    fprintf(fid, 'set ylabel "Energy (eV)"\n');
    fprintf(fid, 'set ylabel offset 1.5,0\n');
    fprintf(fid, 'set xrange [%f:%f]\n', kpoints_l(1), kpoints_l(end));
    fprintf(fid, 'set yrange [%f:%f]\n', E_cut(1), E_cut(2));

    % Set the x-axis tick labels at k-point locations
    fprintf(fid, 'set xtics (');
    for i = 1:length(kpoints_l) - 1
        fprintf(fid, '"%s" %f,', kpoints_name(i), kpoints_l(i));
    end
    fprintf(fid, '"%s" %f)\n', kpoints_name(end), kpoints_l(end));

    % Add arrows to indicate k-point boundaries
    for i = 2:length(kpoints_l) - 1
        fprintf(fid, 'set arrow from %f, %f to %f, %f nohead\n', kpoints_l(i), E_cut(1), kpoints_l(i), E_cut(2));
    end

    % Plot the data
    fprintf(fid, '# plot ''BAND.dat'' u 1:2 w lp lw 2 pt 7 ps 0.2\n');
    fprintf(fid, 'splot "%s.dat" u 1:2:19 w lp lw 2 pt 7 ps 0.2 palette\n', title);
    
    % Set color range for the plot
    fprintf(fid, 'set cbrange [-1:1]\n');
    fprintf(fid, 'unset colorbox\n');

    % Close the file
    fclose(fid);
    
    disp('gnuplot script generated successfully.');
end
