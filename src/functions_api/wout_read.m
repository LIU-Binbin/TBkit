function [orbL, elementL, quantumL] = wout_read(filename, POSCAR_name,options)
%WOUT_READ Parse wannier90.wout files to extract orbital and quantum number data
%
%   Inputs:
%       FILENAME      - Path to wannier90.wout file (default: 'wannier90.wout')
%       POSCAR_NAME   - Path to POSCAR structure file (default: 'POSCAR')
%
%   Outputs:
%       orbL          - Cell array containing orbital projection data
%       elementL      - elements of the wannier functions
%       quantumL      - Quantum number information (currently placeholder)
arguments
    filename = 'wannier90.wout';
    POSCAR_name = 'POSCAR';
    options.UsePOSCAR_Cordinates = 1;
end

% Read crystal structure from POSCAR
[~, sites, ~, ~, elements, ~] =  POSCAR_read(POSCAR_name);
if ~isempty(elements)
    elements.Properties.RowNames = elements.atom_symbol;
end
% Extract metadata from wannier90 output
[SpinfulFlag, Nbands, projLine] = read_metadata(filename);

% Read projection data block
projectionLines = read_projection_data(filename, projLine, Nbands);

% Return extracted data (parse_projections implementation required)
projectionData = parse_projection_cell(projectionLines);
%
orbL = projectionData(:,1:3);  % Orbital projections to be parsed
orbL_Atom = [[sites.rc1];[sites.rc2];[sites.rc3]]';
idx = find_closest_rows(orbL, orbL_Atom);
%
if options.UsePOSCAR_Cordinates
    orbL = orbL_Atom(idx,:);
end
quantumL = [];              % Placeholder for quantum numbers
TestL = [orbL,quantumL];
if has_duplicate_rows_tol(TestL)
    SpinfulFlag = 1;
end
if SpinfulFlag
    spin = 1/2;
    for i = 1:Nbands
        idxIon =  idx(i);
        elementname = remove_atom_numbers(sites(idxIon).name);
        element_data = table2array(elements(...
            char(elementname), {'atom_number','n'}));
        l = projectionData(i,4);
        mr = projectionData(i,5);
        elementL(i,:) = element_data(1) ;
        tmpquantumL = [element_data(2),projectionData(i,4),lmr2m(l,mr),spin];
        if spin == 1/2
            spin = -1/2;
        else
            spin = 1/2;
        end
        quantumL = [quantumL;tmpquantumL];
    end
else
    for i = 1:Nbands
        idxIon =  idx(i);
        elementname = remove_atom_numbers(sites(idxIon).name);
        element_data = table2array(elements(...
            char(elementname), {'atom_number','n'}));
        elementL(i,:) = element_data(1) ;
        l = projectionData(i,4);
        mr = projectionData(i,5);
        tmpquantumL = [element_data(2),projectionData(i,4),lmr2m(l,mr),nan];
        quantumL = [quantumL;tmpquantumL];
    end
end
%


end

function hasDuplicate = has_duplicate_rows_tol(matrix, tol)
if nargin < 2
    tol = 1e-10; % 默认容差
end

sorted_matrix = sortrows(matrix);
row_diffs = diff(sorted_matrix, 1, 1);
hasDuplicate = any(all(abs(row_diffs) < tol, 2));
end

function m = lmr2m(l,mr)
switch l
    case 0
        m = 0;
    case 1
        switch mr
            case 1
                m = 0;
            case 2
                m = 1;
            case 3
                m = -1;
        end
    case 2
        switch mr
            case 1
                m = 0;
            case 2
                m = 1;
            case 3
                m = -1;
            case 4
                m = 2;
            case 5
                m = -2;
        end
    case 3
        switch mr
            case 1
                m = 0;
            case 2
                m = 1;
            case 3
                m = -1;
            case 4
                m = 2;
            case 5
                m = -2;
            case 6
                m = 3;
            case 7
                m = -3;
        end
    case -1
    case -2
    case -3
    case -4
    case -5
end
end
function idx = find_closest_rows(orbL, orbL_atom)
%FIND_CLOSEST_ROWS Find closest row indices between two matrices
%   IDX = FIND_CLOSEST_ROWS(ORBL, ORBL_ATOM) returns the row indices in
%   ORBL_ATOM that are closest to each row in ORBL based on Euclidean distance
%
%   Inputs:
%       orbL      - M x K matrix of vectors (query set)
%       orbL_atom - N x K matrix of vectors (reference set)
%
%   Output:
%       idx       - M x 1 vector of row indices in orbL_atom
%
%   Algorithm:
%       Computes pairwise Euclidean distances efficiently using matrix
%       operations without loops. For large datasets, this vectorized
%       approach is significantly faster than iterative methods.

% 参数验证
if size(orbL, 2) ~= size(orbL_atom, 2)
    error('两个矩阵必须有相同的列数');
end

% 高效计算欧氏距离平方
% 展开距离公式: ||a-b||^2 = ||a||^2 + ||b||^2 - 2*a*b'
dist_sq = sum(orbL.^2, 2) + sum(orbL_atom.^2, 2)' - 2 * (orbL * orbL_atom');

% 为每个orbL行找到最小距离的索引
[~, idx] = min(dist_sq, [], 2);
end


function [SpinfulFlag, Nbands, projLine] = read_metadata(filename)
%READ_METADATA Extract key parameters from wannier90 output file
%   [SPINFLAG, NBANDS, PROJLINE] = READ_METADATA(FILENAME) locates:
%       - Spinor phase inclusion flag
%       - Number of Wannier functions
%       - Starting line of PROJECTIONS section
%
%   Early termination when all metadata is found (typically within first 100 lines)

% Initialize with placeholder values
SpinfulFlag = [];
Nbands = [];
projLine = [];

fid = fopen(filename, 'r');
if fid == -1
    error('File open failed: %s', filename);
end

try
    lineNum = 0;
    while ~feof(fid)
        line = fgetl(fid);
        lineNum = lineNum + 1;

        % Detect spinor phase inclusion
        if isempty(SpinfulFlag) && contains(line, "Include phase for spinor WFs")
            tokens = split_line(line);
            if any(strcmp(tokens, 'T'))
                SpinfulFlag = true;
            else
                SpinfulFlag = false;
            end

        end

        % Extract number of Wannier functions
        if isempty(Nbands) && contains(line, 'Number of Wannier Functions')
            tokens = split_line(line);
            for i = 1:numel(tokens)
                num = str2double(tokens{i});
                if ~isnan(num)
                    Nbands = num;
                    break;
                end
            end
        end

        % Locate PROJECTIONS section start
        if isempty(projLine) && contains(line, 'PROJECTIONS')
            projLine = lineNum;
        end

        % Early exit when all metadata collected
        if ~isempty(SpinfulFlag) && ~isempty(Nbands) && ~isempty(projLine)
            break;
        end
    end
catch ME
    fclose(fid);
    rethrow(ME);
end
fclose(fid);

% Validate critical metadata
if isempty(Nbands)
    error('Band count not found in %s', filename);
end
if isempty(projLine)
    error('PROJECTIONS section not found in %s', filename);
end
end

function projectionLines = read_projection_data(filename, projLine, Nbands)
%READ_PROJECTION_DATA Extract orbital projection lines from file
%   LINES = READ_PROJECTION_DATA(FILENAME, PROJLINE, NBANDS) reads Nbands lines
%   starting at PROJLINE+6 (skipping PROJECTIONS header)
%
%   Implements efficient line-by-line reading without loading entire file

% Calculate target line range [startLine, endLine]
startLine = projLine + 6;
endLine = startLine + Nbands - 1;

fid = fopen(filename, 'r');
if fid == -1
    error('File open failed: %s', filename);
end

try
    % Preallocate cell array for efficiency
    projectionLines = cell(1, Nbands);
    lineNum = 0;
    linesRead = 0;

    while ~feof(fid) && linesRead < Nbands
        line = fgetl(fid);
        lineNum = lineNum + 1;

        % Skip lines before target range
        if lineNum < startLine
            continue;
        end

        % Capture lines within target range
        if lineNum <= endLine
            linesRead = linesRead + 1;
            projectionLines{linesRead} = line;
        else
            break;
        end
    end

    % Validate line count matches expected bands
    if linesRead ~= Nbands
        error('Read %d projection lines, expected %d', linesRead, Nbands);
    end
catch ME
    fclose(fid);
    rethrow(ME);
end
fclose(fid);

end

function tokens = split_line(line, ~)
%SPLIT_LINE Tokenize string while removing empty elements
%   TOKENS = SPLIT_LINE(LINE) splits input using whitespace and common delimiters,
%   returning non-empty tokens in cell array
%
%   Optimized for parsing wannier90 output formatting

% Define standard delimiters for wannier90 files
delimiters = {' ', '|', ':'};

% Perform split and filter empty elements
tokens = strsplit(line, delimiters);
tokens = tokens(~cellfun('isempty', tokens));
end

function matrix = parse_projection_cell(cell_array)
%PARSE_PROJECTION_CELL Convert formatted string cells to numerical matrix
%   MATRIX = PARSE_PROJECTION_CELL(CELL_ARRAY) processes a cell array of
%   strings formatted as ' | 0.74 0.00 ... 1.0 |' into a numerical matrix
%
%   Input:
%       cell_array - Cell array of formatted strings
%
%   Output:
%       matrix     - Numerical matrix with extracted values
%
%   Processing:
%       1. Uses textscan for efficient numerical extraction
%       2. Handles variable whitespace and delimiter characters
%       3. Validates consistent column counts
%       4. Optimized for large datasets

% Validate input
if ~iscell(cell_array)
    error('Input must be a cell array');
end

% Preallocate results matrix
num_rows = numel(cell_array);
results = cell(num_rows, 1);

% Process each cell element
for i = 1:num_rows
    % Extract numbers ignoring non-numeric characters
    nums = textscan(cell_array{i}, '%f', 'Delimiter', ' |', ...
        'MultipleDelimsAsOne', true);

    % Convert to row vector
    results{i} = nums{1}';
end

% Determine matrix dimensions
num_cols = cellfun('length', results);
if any(num_cols ~= num_cols(1))
    error('Inconsistent number of columns in input strings');
end

% Convert to final matrix
matrix = vertcat(results{:});
end

function cleaned = remove_atom_numbers(atom_labels)
%REMOVE_ATOM_NUMBERS Remove numeric suffixes from chemical element labels
%   CLEANED = REMOVE_ATOM_NUMBERS(ATOM_LABELS) processes a cell array of
%   chemical element labels (e.g., {'Te1', 'O2', 'Si'}) and removes any
%   trailing numeric characters, returning only the elemental symbol.
%
%   Input:
%       atom_labels - Cell array of strings with element labels
%
%   Output:
%       cleaned     - Cell array of cleaned element symbols
%
%   Examples:
%       remove_atom_numbers({'Te1', 'O2', 'Si'}) returns {'Te', 'O', 'Si'}
%       remove_atom_numbers('Ca12') returns 'Ca'
%
%   Algorithm:
%       Uses regular expressions to identify and remove trailing digits
%       while preserving the elemental symbol (1-2 alphabetic characters)

% 处理单个字符串输入
if ischar(atom_labels)
    cleaned = regexprep(atom_labels, '(\d+$)', '');
    return;
end

% 验证输入类型
if ~iscell(atom_labels)
    error('Input must be a cell array or string');
end

% 预分配结果数组
cleaned = cell(size(atom_labels));

% 使用正则表达式去除末尾数字
for i = 1:numel(atom_labels)
    % 匹配并删除元素符号后的数字部分
    cleaned{i} = regexprep(atom_labels{i}, '(\d+$)', '');
end
end