function H_hr = add_empty_one(H_hr, vector,options)
%ADD_EMPTY_ONE Add empty hopping terms to the Hamiltonian structure
%   This function extends the Hamiltonian structure H_hr by adding new empty
%   hopping terms specified by the input vector(s). It handles both matrix
%   and sparse representations while maintaining symmetry constraints.
%
%   Inputs:
%       H_hr    - Hamiltonian structure containing:
%                 vectorL: existing hopping vectors
%                 vectorhopping: flag for vector hopping mode
%                 Type: representation type ('mat', 'sparse', or 'list')
%                 WAN_NUM: number of Wannier orbitals
%                 coe/num: flags for symbolic/numeric coefficients
%                 overlap: flag for overlap matrix inclusion
%       vector  - Nx3 matrix of hopping vectors to add ([dx, dy, dz])
%
%   Output:
%       H_hr    - Updated Hamiltonian structure with new empty entries

% Handle vector hopping mode (block diagonal expansion)
arguments
    H_hr 
    vector double
    options.OverlapTest = true;
end
if H_hr(1).overlap && options.OverlapTest
    H_hr(1) = add_empty_one(H_hr(1),vector,'OverlapTest',false);
    H_hr(2) = add_empty_one(H_hr(2),vector,'OverlapTest',false);
    return
end

if H_hr.vectorhopping
    % Append new vectors to hopping list
    H_hr.vectorL = [H_hr.vectorL; vector];
    vector_str = cellfun(@mat2str, vector, 'UniformOutput', false);
    H_hr.vectorL_map =  [H_hr.vectorL_map ;vector_str];
    % Get number of new vectors
    nvector = size(vector,1);
    
    % Expand block diagonal matrices for A/B components
    H_hr.AvectorL = blkdiag(H_hr.AvectorL, eye(nvector));
    H_hr.BvectorL = blkdiag(H_hr.BvectorL, eye(nvector));
    
    % Special handling for C matrix: split and extend both halves
    H_hr.CvectorL = [blkdiag(H_hr.CvectorL(1:end/2,:), eye(nvector));...
                    blkdiag(H_hr.CvectorL(end/2+1:end,:), eye(nvector))];
    return;
end

% Create Map
if isempty(H_hr.vectorL_map)
    % 将矩阵按行分割成元胞数组
    rows = num2cell(H_hr.vectorL, 2);
    % 对每一行应用 mat2str 函数
    H_hr.vectorL_map  = cellfun(@mat2str, rows, 'UniformOutput', false);
end
vector_str = cellfun(@mat2str, num2cell(vector, 2), 'UniformOutput', false);
[New_vector_str,ia] = setdiff(vector_str,H_hr.vectorL_map);
New_vector = vector(ia,:);
N_New_vector = numel(ia);
H_hr.vectorL_map = [H_hr.vectorL_map;New_vector_str];
NRPTS_new = (H_hr.NRPTS +1) : (H_hr.NRPTS )+N_New_vector;
H_hr.vectorL = [H_hr.vectorL;New_vector];

switch H_hr.Type
    case 'mat'  % Full matrix representation
        if H_hr.coe
            H_hr.HcoeL(:,:,NRPTS_new) = zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,N_New_vector,'sym');
        end
        if H_hr.num
            H_hr.HnumL(:,:,NRPTS_new) = zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,N_New_vector,'double');
        end
    case 'sparse'  % Sparse matrix representation
        if H_hr.coe
            newMat{1} = sym(sparse(H_hr.WAN_NUM));
            H_hr.HcoeL{NRPTS_new} = newMat(ones(1,N_New_vector));
        end
        if H_hr.num
            newMat{1} = (sparse(H_hr.WAN_NUM));
            H_hr.HnumL{NRPTS_new} = newMat(ones(1,N_New_vector));
        end
    case 'list'  % List-based representation
        if H_hr.coe
            H_hr.HcoeL(NRPTS_new,1) = sym(0);
        end
        if H_hr.num
            H_hr.HnumL(NRPTS_new,1) = 0;
        end

end

end

