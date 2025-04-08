function [WEIGHTCAR_struct_cell, Name_list, width] = WEIGHTCAR_read_dir(mode)
    % WEIGHTCAR_read_dir generates weightcar structures for 'PBAND' modes
    % This function scans the directory for specific files and processes them accordingly.

    import TBkit_tool.*; % Import necessary functions from the TBkit_tool toolbox

    directory = dir(); % Get the list of files in the current directory
    count = 0; % Initialize counter for files

    % Initialize Name_list to store the processed filenames (without the 'PDOS_' or 'PBAND_')
    Name_list = {""}; 

    switch mode
        case 'DOS'
            % Loop through directory to process 'PDOS_' files
            for i = 1:length(directory)
                filename = directory(i).name;
                % Check if the file contains 'PDOS_' and '.dat'
                if contains(filename, 'PDOS_') && contains(filename, '.dat')
                    count = count + 1;
                    % Store the file name after removing 'PDOS_' and '.dat'
                    Name_list{count} = strrep(strrep(filename, '.dat', ''), 'PDOS_', '');
                    % Process the file and store the density data
                    [Pdensity(:, :, count), ~, ~] = WEIGHTCAR_read(filename, -1, 'vaspkit-DOS-silence');
                end
            end
            % Remove the initial empty entry in Name_list
            Name_list(1) = [];
            % Output the 3D density data structure
            WEIGHTCAR_struct_cell = Pdensity;
            [~, width] = size(Pdensity(:, :, 1)); % Get the width of the data for later use
            
        case 'PBAND'
            % Loop through directory to process 'PBAND_' files
            for i = 1:length(directory)
                filename = directory(i).name;
                % Check if the file contains 'PBAND_' and '.dat'
                if contains(filename, 'PBAND_') && contains(filename, '.dat')
                    count = count + 1;
                    % Store the file name after removing 'PBAND_' and '.dat'
                    Name_list{count} = strrep(strrep(filename, '.dat', ''), 'PBAND_', '');
                    % Process the file and store the 3D band structure data
                    [WEIGHTCAR_3d_mat{count}, ~, ~] = WEIGHTCAR_gen(filename, -1, 'vaspkit-band-silence');
                end
            end
            % Remove the initial empty entry in Name_list
            Name_list(1) = [];
            % Get the width of the first 3D matrix
            width = size(WEIGHTCAR_3d_mat{1}, 3);

            % Convert each 3D matrix into a structured format
            for i = 1:count
                WEIGHTCAR_struct_cell{i} = temp_3d_mat2struct(WEIGHTCAR_3d_mat{i}, Name_list{i}, width);
            end
    end
end

function WEIGHTCAR_struct = temp_3d_mat2struct(WEIGHTCAR_3d, Name, width)
    % temp_3d_mat2struct converts a 3D matrix (weightcar data) into a structured format
    % where each slice (along the third dimension) is stored as a field in a structure.

    import TBkit_tool.*; % Import necessary functions from the TBkit_tool toolbox

    % Determine the orbital names based on the width of the data
    if width > 10
        num2orbital_name = orbital_maprule_gen(1); % Use detailed orbital mapping if width > 10
    else
        num2orbital_name = orbital_maprule_gen(0); % Use basic orbital mapping if width <= 10
    end

    % Initialize the structure array for the output
    for i = 1:width
        % For each slice in the third dimension of the 3D matrix:
        % Assign the weightcar data slice to the structure
        WEIGHTCAR_struct(i).WEIGHTCAR = WEIGHTCAR_3d(:,:,i);
        % Create a display name by combining the Name and the orbital name
        WEIGHTCAR_struct(i).displayname = sprintf('%s-%s', Name, num2orbital_name{i});
    end
end
