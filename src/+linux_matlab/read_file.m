function [outstring] = read_file(filename)
    % READ_FILE: Reads a file line by line and returns the content as a cell array.
    %   filename - Full path to the file to read.
    %
    % Returns:
    %   outstring - Cell array containing each line from the file.
    
    fid = fopen(filename, 'r');  % Open the file for reading.
    
    if fid == -1
        error('Could not open the file: %s', filename);  % Error handling if file cannot be opened.
    end
    
    outstring = {};  % Initialize an empty cell array to store lines.
    
    % Read through the file line by line.
    while ~feof(fid)
        line = fgetl(fid);  % Read one line.
        outstring{end + 1} = line;  % Add the line to the cell array.
    end
    
    fclose(fid);  % Close the file after reading.
end
