function kpath_card_gen()
    % KPOINTS Path Card Generator
    % This function generates two k-point path cards for Wannier90 and wt.in
    % files based on the k-points in the KPOINTS file.
    
    %% Step 1: Read the KPOINTS file
    KPOINTS = fopen('KPOINTS', 'r');  % Open the KPOINTS file for reading
    KPOINTS_information = textscan(KPOINTS, '%s', 'Delimiter', '\n');  % Read all lines of KPOINTS into a cell array
    fclose(KPOINTS);  % Close the KPOINTS file
    
    % Extract number of nodes (usually the 2nd line contains the number of nodes)
    nodes = str2double(KPOINTS_information{1}{2});
    
    %% Step 2: Extract k-point coordinates and labels
    headline = 5;  % Starting line for k-points data
    kpoints = [];  % Initialize k-points array
    kpoints_name_pre = {};  % Initialize k-points labels array
    
    % Loop through the k-points data starting from 'headline'
    for i = headline:length(KPOINTS_information{1})
        kpoints_line = strtrim(KPOINTS_information{1}{i});  % Trim whitespace
        kpoints_line = regexp(kpoints_line, '\s+', 'split');  % Split by whitespace
        
        if length(kpoints_line) > 1  % Ensure line has valid data
            kpoints = [kpoints; str2double(kpoints_line(1:3))];  % Append k-point coordinates
            kpoints_name_pre = [kpoints_name_pre; string(kpoints_line{4})];  % Append k-point label
        end
    end
    
    %% Step 3: Generate the Wannier90 k-point path card
    fid = fopen('wannier90kpath_card', 'w');  % Open file for writing
    fprintf(fid, 'begin kpoint_path\n');  % Write header
    
    % Loop over k-points to write them in the desired format
    for i = 1:2:length(kpoints)
        % Write each pair of k-points with their corresponding labels
        fprintf(fid, '%1s %9.5f %9.5f %9.5f %1s %9.5f %9.5f %9.5f\n', ...
                strrep(kpoints_name_pre{i}, 'GAMMA', 'G'), ...
                kpoints(i, :), ...
                strrep(kpoints_name_pre{i+1}, 'GAMMA', 'G'), ...
                kpoints(i+1, :));
    end
    
    fprintf(fid, 'end kpoint_path\n');  % Write footer
    fclose(fid);  % Close the file
    
    %% Step 4: Generate the wt.in k-point path card
    fid = fopen('wt_in_kpath_card', 'w');  % Open file for writing
    fprintf(fid, 'KPATH_BULK            ! k point path\n');  % Write header
    fprintf(fid, '%d\n', length(kpoints) / 2);  % Write number of k-point pairs
    
    % Loop over k-points to write them in the wt.in format
    for i = 1:2:length(kpoints)
        % Write each pair of k-points with their corresponding labels
        fprintf(fid, '%1s %9.5f %9.5f %9.5f %1s %9.5f %9.5f %9.5f\n', ...
                strrep(kpoints_name_pre{i}, 'GAMMA', 'G'), ...
                kpoints(i, :), ...
                strrep(kpoints_name_pre{i+1}, 'GAMMA', 'G'), ...
                kpoints(i+1, :));
    end
    
    fclose(fid);  % Close the file
end
