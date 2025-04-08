function [outstring, outnumber] = grep(filename, pattern, mode)
    % GREP: Simulates the Unix 'grep' command to search for a pattern in a file.
    %   filename  - Full path to the file to be read.
    %   pattern   - The string or regular expression to match in the file.
    %   mode      - 'aloud' (default) to print the results, or 'silent' to only return them.
    %
    % Returns:
    %   outstring - Cell array of matched lines from the file.
    %   outnumber - Array of line numbers where the pattern was found.

    if nargin < 3
        mode = 'aloud';  % Default mode is 'aloud' to print matching lines.
    end
    
    fid = fopen(filename, 'r');  % Open the file for reading.
    
    if fid == -1
        error('Could not open the file: %s', filename);  % Error handling if file cannot be opened.
    end
    
    line_number = 0;  % Initialize line counter.
    outstring = {};    % Initialize an empty cell array for matched lines.
    outnumber = [];    % Initialize an empty array for line numbers.

    % Read through the file line by line.
    while ~feof(fid)
        line = fgetl(fid);  % Read one line.
        line_number = line_number + 1;  % Increment the line number.
        
        % Search for the pattern within the line.
        matched = strfind(line, pattern);  
        
        % If a match is found, store the line and line number.
        if ~isempty(matched)
            if strcmp(mode, 'aloud')
                fprintf('%d: %s\n', line_number, line);  % Print matching lines if mode is 'aloud'.
            end
            outstring{end + 1} = line;   % Add the line to the output cell array.
            outnumber = [outnumber; line_number];  % Store the line number of the match.
        end
    end
    
    fclose(fid);  % Close the file after reading.
end
