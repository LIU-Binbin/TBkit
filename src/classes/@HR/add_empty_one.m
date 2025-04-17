function H_hr = add_empty_one(H_hr, vector)
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
if H_hr.vectorhopping
    % Append new vectors to hopping list
    H_hr.vectorL = [H_hr.vectorL; vector];
    
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

% Regular mode: process vectors individually
for i = 1:size(vector,1)
    vector_single = vector(i,:);
    
    % Check for existing vectors (skip duplicates)
    try
        % Duplicate check logic:
        % - For non-overlap case: check main vector list
        % - For overlap case: check both main and overlap lists
        if (ismember(vector_single,H_hr.vectorL,'rows') && ~H_hr.overlap) || ...
                (ismember(vector_single,H_hr.vectorL,'rows') ...
                && ismember(vector_single,H_hr.vectorL_overlap,'rows') && H_hr.overlap)
            continue;
        end
    catch
        % Gracefully handle potential dimension mismatches
    end
    
    % Update number of hopping terms
    NRPTS_new = H_hr.NRPTS + 1;
    %H_hr.NRPTS = NRPTS_new;
    
    % Add new vector to list
    H_hr.vectorL(NRPTS_new,:) = vector_single;
    
    % Initialize matrix elements based on representation type
    switch H_hr.Type
        case 'mat'  % Full matrix representation
            if H_hr.coe
                H_hr.HcoeL(:,:,NRPTS_new) = sym(zeros(H_hr.WAN_NUM));
            end
            if H_hr.num
                H_hr.HnumL(:,:,NRPTS_new) = zeros(H_hr.WAN_NUM);
            end
            
            % Handle overlap matrices if required
            if H_hr.overlap
                H_hr.vectorL_overlap(NRPTS_new,:) = vector_single;
                if H_hr.coe
                    H_hr.ScoeL(:,:,NRPTS_new) = sym(zeros(H_hr.WAN_NUM));
                end
                if H_hr.num
                    H_hr.SnumL(:,:,NRPTS_new) = zeros(H_hr.WAN_NUM);
                end
            end
            
        case 'sparse'  % Sparse matrix representation
            if H_hr.coe
                H_hr.HcoeL{NRPTS_new} = sym(sparse(H_hr.WAN_NUM));
            end
            if H_hr.num
                H_hr.HnumL{NRPTS_new} = sparse(H_hr.WAN_NUM);
            end
            
            if H_hr.overlap
                H_hr.vectorL_overlap(NRPTS_new,:) = vector_single;
                if H_hr.coe
                    H_hr.ScoeL{NRPTS_new} = sym(sparse(H_hr.WAN_NUM));
                end
                if H_hr.num
                    H_hr.SnumL{NRPTS_new} = sparse(H_hr.WAN_NUM);
                end
            end
            
        case 'list'  % List-based representation
            if H_hr.coe
                H_hr.HcoeL(NRPTS_new,1) = sym(0);
            end
            if H_hr.num
                H_hr.HnumL(NRPTS_new,1) = 0;
            end
            
            if H_hr.overlap
                H_hr.ScoeL(NRPTS_new,1) = sym(0);
                H_hr.SnumL(NRPTS_new,1) = 0;
            end
    end
end
end

