function write_file(string_list, filename)
    % WRITE_FILE: Writes the contents of string_list to the specified file.
    %   string_list - A cell array of strings to be written to the file.
    %   filename - The name or path of the file to write to.

    % Check if filename is a string
    if ~ischar(filename) && ~isstring(filename)
        error('Filename must be a string or character array');
    end
    
    % Check if string_list is a cell array of strings
    if ~iscell(string_list) || any(~cellfun(@ischar, string_list))
        error('string_list must be a cell array of strings');
    end
    
    % Open file for writing, and handle errors if file cannot be opened
    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not open file %s for writing', filename);
    end
    
    % Write each line to the file
    for i = 1:length(string_list)
        fprintf(fid, '%s\n', string_list{i});
    end
    
    % Close the file
    fclose(fid);
end
