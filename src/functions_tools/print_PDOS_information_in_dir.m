function print_PDOS_information_in_dir(Pdos_namelist, width)
    % PRINT_PDOS_INFORMATION_IN_DIR Displays pDOS data information.
    %
    %   print_PDOS_INFORMATION_IN_DIR(Pdos_namelist, width)
    %
    % Inputs:
    %   Pdos_namelist - Cell array of pDOS data filenames.
    %   width          - Integer indicating the number of columns.
    %
    % Outputs:
    %   None
    %
    % Note:
    %   Ensure that Pdos_namelist is a cell array of strings and width is an integer.
    %
    % Author: [Your Name]
    % Date: [Date]
    % Version: 1.0
    % License: [License Information]
    % Contact: [Your Contact Information]

    % Validate input parameters
    if nargin < 2
        error('Insufficient input arguments. Please provide Pdos_namelist and width.');
    end

    % Display the list of pDOS data files
    fprintf('We have these pDOS data files; the mapping is:\n');
    for i = 1:length(Pdos_namelist)
        fprintf(' %d  :  %s\n', i, Pdos_namelist{i});
    end

    % Display the format based on the width
    fprintf('And the format is: ');
    if width > 10
        fprintf('s py pz px dxy dyz dz2 dxz dx2-y2 fy3x2 fxyz fyz2 fz3 fxz2 fzx2 fx3 tot\n');
        fprintf('1  2  3  4   5   6   7   8      9    10    11    12  13   14   15   16  17\n');
    else
        fprintf('s py pz px dxy dyz dz2 dxz dx2-y2 tot\n');
        fprintf('1  2  3  4   5   6   7   8      9  10\n');
    end
end
