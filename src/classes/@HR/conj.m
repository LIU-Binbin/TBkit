function H_hr = conj(H_hr)
% CONJ Complex conjugate of HR object in momentum space representation
%
%   H_conj = CONJ(H_hr) computes the complex conjugate Hamiltonian that
%   satisfies the time-reversal symmetry condition in momentum space:
%       H_conj(k) = H_hr(-k)^*
%
%   The real-space implementation:
%   1. Dualizes the R-space vectors (r -> -r) using dualize()
%   2. Applies complex conjugation to all coupling matrices
%   3. Preserves original R-vector ordering through index mapping
%
%   Input:
%       H_hr - HR object with Hamiltonian data
%           Required properties:
%               dualize() method implementation
%               Duality_vector_dist permutation vector
%               HnumL/HcoeL: Hamiltonian matrices
%               Type: Storage format ('mat' or 'sparse')
%
%   Output:
%       H_hr - New HR object with conjugated Hamiltonian
%
%   Example:
%       H_conj = conj(H_hr); % Get time-reversed Hamiltonian
%
%   See also DUALIZE, CTRANSPOSE, AUTOHERMI

    % Validate object state
    try
        H_hr = H_hr.dualize();
    catch
        error('Dualize method failed. Verify HR object implementation.');
    end
    
    % Check permutation vector existence
    if ~isprop(H_hr, 'Duality_vector_dist') || isempty(H_hr.Duality_vector_dist)
        error('Missing Duality_vector_dist. Verify dualize() implementation.');
    end
    
    % Validate permutation indices
    if any(H_hr.Duality_vector_dist > H_hr.NRPTS) || any(H_hr.Duality_vector_dist < 1)
        error('Invalid Duality_vector_dist indices. Check dualize() mapping.');
    end
    
    % Matrix conjugation operation
    if strcmp(H_hr.Type, 'mat')
        % Dense matrix conjugation
        if  H_hr.coe
            H_hr.HcoeL = conj(H_hr.HcoeL(:,:,H_hr.Duality_vector_dist));
        end
        if H_hr.num
            H_hr.HnumL = conj(H_hr.HnumL(:,:,H_hr.Duality_vector_dist));
        end
    
    elseif strcmp(H_hr.Type, 'sparse')
        % Sparse matrix conjugation
        if H_hr.coe
            H_hr.HcoeL = cellfun(@conj, H_hr.HcoeL(H_hr.Duality_vector_dist), 'UniformOutput', false);
        end
        if H_hr.num
            H_hr.HnumL = cellfun(@conj, H_hr.HnumL(H_hr.Duality_vector_dist), 'UniformOutput', false);
        end
    elseif strcmp(H_hr.Type, 'list')
        % add list !!!
        % Process list-format Hamiltonian
        NRPTS_ = H_hr.NRPTS; % Get original number of lattice vectors
        vectorList = H_hr.vectorL; % Extract original vector list
        % Create opposite vectors by negating spatial components
        vectorList_oppo(:,1:H_hr.Dim) = -vectorList(:,1:H_hr.Dim);
    
        % Handle special 5-column format (3 spatial + 2 extra dimensions)
        if size(vectorList,2) == 5
            % Swap last two columns for opposite vectors
            vectorList_oppo(:,H_hr.Dim+1) = vectorList(:,H_hr.Dim+1);
            vectorList_oppo(:,H_hr.Dim+2) = vectorList(:,H_hr.Dim+2);
        end
    
        % Initialize duality mapping vector
        Duality_vector_dist = zeros(NRPTS_,1);
    
        % Process each lattice vector to find/create its dual
        for i = 1:NRPTS_
            % Get current dual candidate vector
            vector_tmp_oppo = vectorList_oppo(i,:);
    
            % Find existing dual vector in original list
            [~,j] = ismember(vector_tmp_oppo, H_hr.vectorL, 'rows');
    
            % Add new empty entry if dual vector not found
            if j == 0
                H_hr = H_hr.add_empty_one(vector_tmp_oppo);  % Expand vector list
                j = H_hr.NRPTS;                             % Get new index
                Duality_vector_dist(j) = i;            % Record reverse mapping
            end
    
            % Store duality mapping: original -> dual
            Duality_vector_dist(i) = j;
        end
        if H_hr.num
            H_hr.HnumL = (H_hr.HnumL(H_hr.Duality_vector_dist));
        end
        if H_hr.coe
            H_hr.HcoeL = (H_hr.HcoeL(H_hr.Duality_vector_dist));
        end
    
    end
end