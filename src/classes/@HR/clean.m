function H_hr = clean(H_hr, WANNUM)
% CLEAN Reset Hamiltonian and overlap matrices while preserving structure
%
%   H_hr = CLEAN(H_hr) resets all Hamiltonian (and overlap) matrices to zero 
%   while maintaining the original Wannier orbital count and structure
%
%   H_hr = CLEAN(H_hr, WANNUM) resizes matrices to specified orbital count
%
%   Inputs:
%       H_hr    - HR object to be cleaned
%       WANNUM  - (Optional) New number of Wannier orbitals
%
%   Output:
%       H_hr    - Cleaned HR object with zero-initialized matrices
%
%   Features:
%       - Preserves R-vector list and symmetry information
%       - Handles both numeric and symbolic coefficients
%       - Maintains overlap matrices if present
%       - Supports matrix and sparse storage formats
%
%   See also HR, AUTOHERMI, CHARALIZE

    arguments
        H_hr HR
        WANNUM {mustBeInteger,mustBePositive} = H_hr.WAN_NUM
    end

    % Update Wannier orbital count
    H_hr.WAN_NUM = WANNUM;
    
    % Handle different storage formats
    switch H_hr.Type
        case 'mat'
            % Dense matrix format
            reset_dense_matrices();
            
        case 'Sparse'
            % Sparse matrix format
            reset_sparse_matrices();
            
        otherwise
            error('Unsupported matrix type: %s', H_hr.Type);
    end

    % Nested helper functions
    function reset_dense_matrices()
        % Reset Hamiltonian matrices
        H_hr.HnumL = zeros(WANNUM, WANNUM, H_hr.NRPTS);
        
        % Handle symbolic coefficients
        if H_hr.coe
            H_hr.HcoeL = sym(H_hr.HnumL);
        else
            H_hr.HcoeL = sym([]);
        end
        
        % Process overlap matrices if present
        if isfield(H_hr, 'overlap') && H_hr.overlap
            H_hr.SnumL = zeros(WANNUM, WANNUM, H_hr.NRPTS);
            
            if H_hr.coe
                H_hr.ScoeL = sym(H_hr.SnumL);
            else
                H_hr.ScoeL = sym([]);
            end
        end
    end

    function reset_sparse_matrices()
        % Reset sparse storage
        warning('Sparse matrix cleaning not fully implemented');
        
        % Initialize sparse structure
        H_hr.HnumL = cell(1, H_hr.NRPTS);
        if H_hr.coe
            H_hr.HcoeL = cell(1, H_hr.NRPTS);
        end
        
        % Create empty sparse matrices for each R-point
        for i = 1:H_hr.NRPTS
            H_hr.HnumL{i} = sparse(WANNUM, WANNUM);
            if H_hr.coe
                H_hr.HcoeL{i} = sym(sparse(WANNUM, WANNUM));
            end
        end
        
        % Handle overlap matrices
        if isfield(H_hr, 'overlap') && H_hr.overlap
            H_hr.SnumL = cell(1, H_hr.NRPTS);
            if H_hr.coe
                H_hr.ScoeL = cell(1, H_hr.NRPTS);
            end
        end
    end
end