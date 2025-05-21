function [dataArray, NRPT_list, NRPTS, NUM_WAN] = hrdat_read(filename)
%HRDAT_READ Reads Wannier90 Hamiltonian data from _hr.dat file
%   [DATAARRAY, NRPT_LIST, NRPTS, NUM_WAN] = HRDAT_READ(FILENAME) reads:
%       - NUM_WAN: Number of Wannier functions
%       - NRPTS: Number of real-space points
%       - NRPT_LIST: Degeneracy weights of k-points
%       - DATAARRAY: Hopping parameters [i, j, nx, ny, nz, real_part, imag_part]
%
%   If no filename specified, defaults to 'wannier90_hr.dat'

    %% Handle default filename
    if nargin < 1
        filename = 'wannier90_hr.dat';
    end

    %% Read header information (NUM_WAN and NRPTS)
    fileID = fopen(filename, 'r');
    if fileID == -1
        error('Failed to open file: %s', filename);
    end
    
    % Read NUM_WAN (number of Wannier functions) from line 2
    % Read NRPTS (number of real-space points) from line 3
    headerFormat = '%d%*s';  % Read integer, ignore remaining content
    headerData = textscan(fileID, headerFormat, 2,...
                         'HeaderLines', 1,...
                         'Delimiter', ' ',...
                         'MultipleDelimsAsOne', true);
    
    NUM_WAN = headerData{1}(1);
    NRPTS = double(headerData{1}(2));
    fclose(fileID);

    %% Read NRPT_LIST (k-point weights)
    fileID = fopen(filename, 'r');
    
    % Calculate number of full lines (15 weights per line)
    weightsPerLine = 15;
    fullLines = floor(NRPTS/weightsPerLine);
    remainingWeights = mod(NRPTS, weightsPerLine);
    
    % Read weight data blocks
    weightFormat = repmat('%f', 1, weightsPerLine);
    weightData = textscan(fileID, weightFormat, fullLines,...
                         'HeaderLines', 3,...
                         'Delimiter', ' ',...
                         'MultipleDelimsAsOne', true);
    
    % Read remaining weights if exists
    if remainingWeights > 0
        partialFormat = repmat('%f', 1, remainingWeights);
        partialData = textscan(fileID, partialFormat, 1,...
                              'Delimiter', ' ',...
                              'MultipleDelimsAsOne', true);
        weightData = [weightData, partialData];
    end
    
    % Flatten cell array to column vector
    NRPT_list = cell2mat(weightData(:));
    fclose(fileID);

    %% Read hopping parameters
    fileID = fopen(filename, 'r');
    
    % Skip header and weight sections (4 + fullLines + ceil(remainingWeights/15))
    skipLines = 3 + double(fullLines) + (remainingWeights > 0);
    hoppingFormat = '%f %f %f %f %f %f %f';  % [i, j, nx, ny, nz, real, imag]
    
    dataArray = textscan(fileID, hoppingFormat,...
                        'HeaderLines', skipLines,...
                        'Delimiter', ' ',...
                        'MultipleDelimsAsOne', true);
    
    fclose(fileID);
end