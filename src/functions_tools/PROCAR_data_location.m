function line = PROCAR_data_location(ion, kpoint, band, NUMFLAG, mode, ionsnum, bandsnum, SOC_flag)
    % PROCAR_data_location - Calculates the line number of a specific data entry in the PROCAR file.
    %
    % Inputs:
    %   ion        - Ion index for which the data is being retrieved.
    %   kpoint     - K-point index in the Brillouin zone.
    %   band       - Band index.
    %   NUMFLAG    - Data component selector (0 for TOT, 1 for MX, etc.).
    %   mode       - Access mode ('f' for standard access, can be extended).
    %   ionsnum    - Total number of ions in the system.
    %   bandsnum   - Total number of bands in the system.
    %   SOC_flag   - Flag for spin-orbit coupling (default: 1).
    %
    % Output:
    %   line       - Line number where the requested data can be found in PROCAR.
    
    % Default arguments if not provided
    if nargin < 8
        SOC_flag = 1;  % Default: SOC flag is 1 (on)
    end
    if nargin < 6
        % Read PROCAR header information if not provided
        filename = 'PROCAR_information';
        delimiter = ' ';
        startRow = 2;
        formatSpec = '%*s%*s%*s%f%*s%*s%*s%f%*s%*s%*s%f%[^\n\r]';
        
        % Try to open and read the file
        try
            fileID = fopen(filename, 'r');
            textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false);
            fclose(fileID);
            
            PROCARinformatiobn = [dataArray{1:end-1}];
            
            % Extract header information
            kpointsnum = PROCARinformatiobn(1);
            bandsnum = PROCARinformatiobn(2);
            ionsnum = PROCARinformatiobn(3);
        catch
            error('Could not read PROCAR file or open PROCAR_information.');
        end
        
        clearvars filename delimiter startRow formatSpec fileID dataArray ans;
    end
    
    if nargin < 5
        mode = 'f';  % Default mode
    end
    
    %% Calculate data line number based on SOC flag
    if strcmp(mode, 'f')
        if SOC_flag == 1
            % Spin-orbit coupling is included
            kpoints_round = bandsnum * ionsnum * 4;
            bands_round = ionsnum * 4;
            line = (kpoint - 1) * kpoints_round + (band - 1) * bands_round + NUMFLAG * ionsnum + ion;
        else
            % No spin-orbit coupling
            kpoints_round = bandsnum * ionsnum;
            bands_round = ionsnum;
            line = (kpoint - 1) * kpoints_round + (band - 1) * bands_round + NUMFLAG * ionsnum + ion;
        end
    else
        error('Mode %s not supported yet.', mode);
    end
end
