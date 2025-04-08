function hspiceACTR(filename)
    % HSPICEACTR: Process an HSPICE file by creating a backup, finding lines, 
    % and running PowerShell commands to extract and manipulate content.
    %
    % Input:
    %   filename - The name of the HSPICE file to process.
    
    % Check if the file exists
    if ~isfile(filename)
        error('The specified file does not exist: %s', filename);
    end
    
    % Create a backup of the file
    backup_filename = strcat(filename, '.bk');
    try
        copyfile(filename, backup_filename);
    catch
        error('Failed to create a backup file: %s', backup_filename);
    end
    
    % Find lines containing the specific pattern
    Selectline = win_matlab.FindStrLines(filename, "$&%#");
    
    % Validate the Selectline value
    if isempty(Selectline)
        error('No matching lines found with the pattern "$&%#".');
    end
    
    % Run PowerShell to extract the first part of the file
    PowerShellCmd = sprintf('powershell -command "type %s | Select -First %s > %s.info"', ...
        filename, string(Selectline), filename);
    try
        system(PowerShellCmd);
    catch
        error('PowerShell command failed to execute: %s', PowerShellCmd);
    end
    
    % Run PowerShell to skip specific lines and save the remaining content
    PowerShellCmd = sprintf('powershell -command "Get-Content %s | Select-Object -Skip %s | Set-Content %s.data"', ...
        filename, string(Selectline), filename);
    try
        system(PowerShellCmd);
    catch
        error('PowerShell command failed to execute: %s', PowerShellCmd);
    end
    
    fprintf('Processing complete. Backup and extracted data saved.\n');
end
