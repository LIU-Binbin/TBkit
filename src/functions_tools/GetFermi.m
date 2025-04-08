function Efermi = GetFermi(mode, filename)
% GetFermi - A function to extract the Fermi energy from different types of input files.
%
% Usage:
%   Efermi = GetFermi(mode, filename)
%
% Inputs:
%   mode      - A string that determines the format of the file (e.g., 'vasp', 'qe').
%   filename  - The name of the file to read the Fermi energy from (default is 'DOSCAR').
%
% Outputs:
%   Efermi    - The extracted Fermi energy value.

    % Set the default filename if not provided
    if nargin < 2
        filename = 'DOSCAR';  % Default filename for VASP
    end
    
    % Check the mode for different formats
    if strcmp(mode, 'vasp')
        % VASP format: Read the Fermi energy from the DOSCAR file
        
        % Define parameters for reading the file
        delimiter = ' ';
        startRow = 6;    % Start reading from line 6
        endRow = 6;      % Only read 1 line (line 6)
        formatSpec = '%*s%*s%*s%f%*s%[^\n\r]';  % Format to extract the Fermi energy (4th column)
        
        % Open the file and extract the data
        fileID = fopen(filename, 'r');
        dataArray = textscan(fileID, formatSpec, endRow-startRow+1, ...
            'Delimiter', delimiter, 'MultipleDelimsAsOne', true, ...
            'TextType', 'string', 'HeaderLines', startRow-1, ...
            'ReturnOnError', false, 'EndOfLine', '\r\n');
        fclose(fileID);
        
        % Store the Fermi energy value (extracted from the 4th column)
        Efermi = dataArray{1};
        
        % Clear variables used for file reading
        clearvars filename delimiter startRow endRow formatSpec fileID dataArray;
        
    elseif strcmp(mode, 'qe')
        % Quantum Espresso format: Extract Fermi energy from the 'scf.out' file
        
        % Search for the line containing 'Fermi' in the 'scf.out' file
        [Efermi_string,~] = grep('scf.out', 'Fermi');
        
        % Extract numerical value of Fermi energy from the string
        Efermi = str2double(regexp(Efermi_string, '\d*\.?\d*', 'match', 'once'));
        
        % Display the Fermi energy for debugging purposes
        disp(['Fermi Energy (QE): ', num2str(Efermi)]);
        
    else
        % Return 0 for unsupported modes
        Efermi = 0;
        disp('Unsupported mode or file format.');
    end
end
