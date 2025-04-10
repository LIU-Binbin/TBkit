function H_htrig = setup(H_htrig, Var_cell, k_cell, mat_cell, silence)
%SETUP Configure an Htrig object with symbolic or numeric hopping terms.
%
%   H_htrig = SETUP(H_htrig, Var_cell, k_cell, mat_cell)
%   Initializes the Htrig object's hopping structure by assigning
%   hopping matrices (mat_cell) to symbolic expressions (Var_cell)
%   associated with momentum symbols (k_cell).
%
%   H_htrig = SETUP(H_htrig, Var_cell, k_cell, mat_cell, silence)
%   Suppresses progress bar output if 'silence' is true.
%
%   Inputs:
%       H_htrig   - An instance of Htrig class
%       Var_cell  - Cell array of symbolic or numeric prefactors
%       k_cell    - Cell array of corresponding k-symbols (momentum-dependent)
%       mat_cell  - Cell array of hopping matrices (numeric or symbolic)
%       silence   - (Optional) If true, disables progress bar (default: false)
%
%   Output:
%       H_htrig   - Updated Htrig object with configured hopping terms
%
%   Example:
%       Var_cell = {cos(kx), sin(kx)};
%       k_cell = {exp(1i*kx), exp(-1i*kx)};
%       mat_cell = {t*eye(2), t*eye(2)};
%       H = setup(H, Var_cell, k_cell, mat_cell);
%
%   Notes:
%       - Htrig.HsymL_trig and Htrig.HcoeL/HnumL will be updated accordingly.
%       - Automatically distinguishes symbolic and numeric modes.
%       - Relies on k_symbol2Kind to manage k-symbol indexing.
%
%   See also: k_symbol2Kind, TBkit_tool_outer.CmdLineProgressBar

if nargin < 5
    silence = false;
end

BASIS_NUM = H_htrig.Basis_num;

if length(Var_cell) ~= length(k_cell) || length(k_cell) ~= length(mat_cell)
    error('Mismatch in lengths of input cell arrays.');
end

if ~silence
    pb = TBkit_tool_outer.CmdLineProgressBar('Setting ');
end

nVar_cell = length(Var_cell);

for i = 1:nVar_cell
    Var = Var_cell{i};
    k_symbol = k_cell{i};
    matcell = mat_cell{i};

    Kind = H_htrig.k_symbol2Kind(k_symbol);
    if isempty(Kind)
        Kind = H_htrig.Kinds + 1;
        H_htrig.HsymL_trig(Kind) = k_symbol;
        H_htrig.HcoeL(:,:,Kind) = sym(zeros(BASIS_NUM, BASIS_NUM));
        H_htrig.HnumL(:,:,Kind) = zeros(BASIS_NUM, BASIS_NUM);
    end

    if ~silence
        pb.print(i, nVar_cell, 'Htrig ...');
    end

    switch class(Var)
        case 'sym'
            H_htrig.HcoeL(:,:,Kind) = H_htrig.HcoeL(:,:,Kind) + matcell * Var;
        case 'double'
            H_htrig.HnumL(:,:,Kind) = H_htrig.HnumL(:,:,Kind) + matcell * Var;
        case 'string'
            % Placeholder: add custom behavior here if needed
        otherwise
            error('Unsupported variable type: %s', class(Var));
    end
end

if ~silence
    pb.delete();
end
end

