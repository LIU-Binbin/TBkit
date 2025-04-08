function Line = TotalLines(filename)
    % TotalLines: Counts the total number of lines in a file.
    % On Windows, it uses the system 'find' command.
    % On Unix-based systems, it uses 'wc -l'.
    %
    % Input:
    %   filename - The path to the file for which to count the lines.
    %
    % Output:
    %   Line - The total number of lines in the file.
    % Example：
    %   filename = 'example.txt';
    %   totalLines = TotalLines(filename);
    %   fprintf('The total number of lines in the file is: %d\n', totalLines);

    % Check if the file exists
    if ~isfile(filename)
        error('The specified file does not exist: %s', filename);
    end
    
    % Initialize Line variable
    Line = 0;

    if ispc  % Windows system
        % Use find command to count lines (ignores empty lines)
        [status, output] = system(['find /c /v "" ', filename]);
        if status ~= 0
            error('Failed to execute system command on Windows.');
        end
        % Parse output: extract number of lines from the find command output
        tokens = strsplit(output, ':');
        Line = str2double(strtrim(tokens{2}));
        
    elseif isunix  % Unix/Linux/Mac system
        % Use wc -l to count the total number of lines in the file
        [status, output] = system(['wc -l < ', filename]);
        if status ~= 0
            error('Failed to execute system command on Unix-based system.');
        end
        % wc -l returns the line count as output, directly convert it to number
        Line = str2double(strtrim(output));
        
    else
        error('Unsupported operating system.');
    end
end
