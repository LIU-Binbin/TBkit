function [WEIGHTCAR, EIGENCAR, KPATH] = WEIGHTCAR_read(filename, seq_list, mode)
% WEIGHTCAR_READ Reads and processes VASP band structure and DOS projection data
% Input Parameters:
%   filename  - Input data file name
%   seq_list  - Orbital indices to sum (default: all orbitals)
%   mode      - Processing mode ('vaspkit-band', 'vaspkit-dos', etc.)
% Output Parameters:
%   WEIGHTCAR - Weight data matrix/array
%   EIGENCAR  - Eigenvalue data
%   KPATH     - K-point path coordinates (band mode only)

% Parameter initialization
arguments
    filename 
    seq_list  = -1;
    mode = 'vaspkit-band';
end

% Initialize output variables
[WEIGHTCAR, EIGENCAR, KPATH] = deal([]);

try
    % Main processing logic
    if contains(lower(mode), 'band')
        processBandData();
    elseif contains(lower(mode), 'dos')
        processDOSData();
    else
        error('Unsupported processing mode: %s', mode);
    end
catch ME
    handleCleanup();
    rethrow(ME);
end

% Final cleanup
handleCleanup();

% ================== Nested Functions ================== %
    function processBandData()
        % Process band structure data
        silentMode = contains(mode, 'silent');
        
        % Read and preprocess data
        [data, nBands] = readAndFilterData();
        
        % Extract K-path and eigenvalues
        KPATH = data(1:size(data,1)/nBands, 1);
        EIGENCAR = reshape(data(:,2), [], nBands).';
        
        % Process projection data
        projData = data(:, 3:end);
        [nKpts, nOrbitals] = deal(size(projData,1)/nBands, size(projData,2));
        
        % Reshape and handle band ordering
        proj3D = permute(reshape(projData, nKpts, nBands, nOrbitals), [2 1 3]);
        proj3D(2:2:end,:,:) = flip(proj3D(2:2:end,:,:), 2);
        
        % Handle orbital selection
        if ~usingDefaultSeq()
            validateOrbitalIndices(nOrbitals);
            WEIGHTCAR = sum(proj3D(:,:,seq_list), 3);
        else
            WEIGHTCAR = proj3D;
        end
        
        % Output information
        if ~silentMode
            printBandHeader(nOrbitals);
        end
    end

    function processDOSData()
        % Process DOS data
        silentMode = contains(mode, 'silent');
        
        % Read data
        data = readmatrix(filename, 'NumHeaderLines', 1);
        EIGENCAR = data(:,1);
        projData = data(:,2:end);
        [~, nOrbitals] = size(projData);
        
        % Handle orbital selection
        if ~usingDefaultSeq()
            validateOrbitalIndices(nOrbitals);
            WEIGHTCAR = sum(projData(:,seq_list), 2);
        else
            WEIGHTCAR = projData;
        end
        
        % Output information
        if ~silentMode
            printDOSHeader(nOrbitals);
        end
    end

    function [data, nBands] = readAndFilterData()
        % Read and filter band data
        fileContent = fileread(filename);
        lines = splitlines(fileContent);
        
        % Remove comment lines and headers
        validLines = ~contains(lines, {'#', 'Band'});
        dataLines = lines(validLines);
        
        % Convert to numeric matrix
        % 超大数据量优化版
        all_data = strjoin(dataLines, ' ');
        data = sscanf(all_data, '%f', [numel(str2num(dataLines{1})), numel(dataLines)])';
        nBands = sum(contains(lines, 'Band'));
        
        % Validate data integrity
        if mod(size(data,1), nBands) ~= 0
            error('Data row count mismatch with number of bands');
        end
    end

    function flag = usingDefaultSeq()
        % Check for default orbital selection
        flag = isequal(seq_list, -1) || isempty(seq_list);
    end

    function validateOrbitalIndices(maxOrbital)
        % Validate orbital indices
        if any(seq_list > maxOrbital) || any(seq_list < 1)
            error('Orbital indices out of range [1-%d]', maxOrbital);
        end
    end

    function printBandHeader(nOrbitals)
        % Display band structure header
        fprintf('Vaspkit 1.2.0 band format supported\n');
        if nOrbitals >= 17
            fprintf('Orbital components: s py pz px dxy dyz dz2 dxz dx2-y2 fy3x2 fxyz fyz2 fz3 fxz2 fzx2 fx3 tot\n');
        else
            fprintf('Orbital components: s py pz px dxy dyz dz2 dxz dx2-y2 tot\n');
        end
    end

    function printDOSHeader(nOrbitals)
        % Display DOS header
        fprintf('Vaspkit 1.2.0 DOS format supported\n');
        if nOrbitals >= 17
            fprintf('Orbital components: s py pz px dxy dyz dz2 dxz dx2-y2 fy3x2 fxyz fyz2 fz3 fxz2 fzx2 fx3 tot\n');
        else
            fprintf('Orbital components: s py pz px dxy dyz dz2 dxz dx2-y2 tot\n');
        end
    end

    function handleCleanup()
        % Cleanup temporary files
        tempFiles = {'temp_Pband.dat', 'temp_Pband.DAT'};
        for f = tempFiles
            if exist(f{1}, 'file'), delete(f{1}); end
        end
    end
end