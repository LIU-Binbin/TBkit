function H_hr = from_wannier90(filename, Type, options)
%FROM_WANNIER90 Construct HR object from Wannier90 _hr.dat file(s)
%   This function imports tight-binding model data from Wannier90 output files,
%   http://www.wanniertools.com/input.html#dense-format-storage
%   supporting both regular Hamiltonian and overlap matrix import.
%
%   H_hr = FROM_WANNIER90(filename, Type, options) creates an HR object from:
%       - Single _hr.dat file (regular Hamiltonian)
%       - Pair of _hr.dat files (Hamiltonian + overlap matrix)
%
%   Input Parameters:
%       filename : Path(s) to wannier90_hr.dat file(s)
%                   - Single string: Regular Hamiltonian
%                   - Cell array {H_file, S_file}: Hamiltonian + Overlap
%       Type     : Data organization type ('mat' matrix form / 'list' sparse form)
%       options  : Struct with processing parameters:
%           - Accuracy: Numerical threshold for element filtering (default: 1e-6)
%           - overlap : Flag for overlap matrix inclusion (default: false)
%
%   Output:
%       H_hr : HR object containing Hamiltonian (and optionally overlap) data
%               with properties:
%               - HnumL/SnumL: Numerical Hamiltonian/Overlap coefficients
%               - vectorL: Lattice vectors
%               - WAN_NUM: Number of Wannier functions
%
%   Processing Features:
%       - Automatic degeneracy handling using NRPT_list
%       - Numerical filtering based on specified accuracy
%       - Support for both matrix and list storage formats
%       - Symmetry flag initialization (default: false)
%
%   See also HR, hrdat_read, sym

arguments
    filename (1,:) {validateFilename} = 'wannier90_hr.dat'  % File(s) path
    Type char {mustBeMember(Type,{'mat','list'})} = 'mat'   % Data storage format
    options.Accuracy double {mustBePositive} = 1e-6         % Numerical threshold
    options.overlap logical = false                         % Overlap matrix flag
    options.OverlapTest logical= true;
end
% Main processing branch for overlap matrix case
if options.overlap
    options.OverlapTest = false;
    options.overlap = false;
    optionscell = namedargs2cell(options);

    % Read Hamiltonian and overlap files
    H_hr(1) = HR.hrdat_read(filename{1},Type,optionscell{:});
    H_hr(2) = HR.hrdat_read(filename{2},Type,optionscell{:});

% Regular Hamiltonian processing (no overlap)    
else
    % Read single Hamiltonian file
    [dataArray, NRPT_list, NRPTS, NUM_WAN] = HR.hrdat_read(filename);
    
    % Matrix format processing
    if strcmp(Type, 'mat')
        [vectorL, HnumL] = processMatrixData(dataArray, NRPT_list, NUM_WAN, NRPTS);
        HcoeL = sym([]);
        
    % List format processing    
    elseif strcmp(Type, 'list')
        [vectorL, HnumL] = processListData(dataArray, NRPT_list, options.Accuracy);
        HcoeL = sym([]);
    end
    
    % Create basic HR object
    H_hr = HR(NUM_WAN, vectorL,...
        'HnumL', HnumL,...
        'HcoeL', HcoeL,...
        'Type', Type,...
        'sym', false);
    
    % Set numerical properties
    H_hr.num = true;
    H_hr.coe = false;
    H_hr.Basis_num = H_hr.WAN_NUM;
end

%% try to read orbital infomation from wout file
if exist("wannier90.wout","file")
    [orbL, elementL, quantumL] = wout_read("wannier90.wout", "POSCAR");
    H_hr.orbL = orbL;
    H_hr.elementL = elementL;
    H_hr.quantumL = quantumL;
end

%% Nested helper functions
    function [vec, mat] = processMatrixData(data, nrpt, numWan, nrpts)
        %PROCESSMATRIXDATA Convert raw data to matrix format
        % Extracts vector components and reshapes Hamiltonian/overlap matrices
        Vec_Fir = data{:, 1};      % a1 direction
        Vec_Sec = data{:, 2};      % a2 direction
        Vec_Thi = data{:, 3};      % a3 direction
        %Orb_fir = dataArray{:, 4};
        %Orb_sec = dataArray{:, 5};
        h_real = data{:, 6};
        h_imag = data{:, 7};
        
        % Reshape vectors
        V_f = reshape(Vec_Fir, numWan*numWan, nrpts);
        V_s = reshape(Vec_Sec, numWan*numWan, nrpts);
        V_t = reshape(Vec_Thi, numWan*numWan, nrpts);
        vec = [V_f(1,:)', V_s(1,:)', V_t(1,:)'];
        
        % Reshape and combine complex numbers
        mat_real = reshape(h_real, numWan, numWan, nrpts);
        mat_imag = reshape(h_imag, numWan, numWan, nrpts);
        mat = mat_real + 1i*mat_imag;
        
        % Normalize if requested
        allone = all(double(nrpt) == 1);% check nrpt all 1
        if ~allone
            for i = 1:nrpts
                mat(:,:,i) = mat(:,:,i)/nrpt(i); % bug?
            end
        end
    end

    function [vec, mat] = processListData(data, nrpt, accuracy)
        %PROCESSLISTDATA Convert raw data to sparse list format
        % Filters small elements and organizes data in sparse representation
        dataMatrix = cell2mat(data(1:7));
        complexData = dataMatrix(:,6) + 1i*dataMatrix(:,7);
        
        % Apply accuracy threshold
        validIdx = abs(complexData) > accuracy;
        vec = dataMatrix(validIdx, 1:3);
        [~, ~, origLabel] = unique(vec(:,1:3), 'rows');
        
        % Normalize using degeneracy factors
        mat = complexData(validIdx) ./ nrpt(origLabel);
    end
end

%% Validation functions
function validateFilename(fname)
%VALIDFILENAME Ensure correct filename input type
    if ~(ischar(fname) || (iscell(fname) && numel(fname)==2))
        error('Filename must be string or 2-element cell array');
    end
end