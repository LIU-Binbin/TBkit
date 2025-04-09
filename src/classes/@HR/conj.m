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
    if ~isfield(H_hr, 'Duality_vector_dist') || isempty(H_hr.Duality_vector_dist)
        error('Missing Duality_vector_dist. Verify dualize() implementation.');
    end
    
    % Validate permutation indices
    if any(H_hr.Duality_vector_dist > H_hr.NRPTS) || any(H_hr.Duality_vector_dist < 1)
        error('Invalid Duality_vector_dist indices. Check dualize() mapping.');
    end

    % Matrix conjugation operation
    if strcmp(H_hr.Type, 'mat')
        % Dense matrix conjugation
        if ~isempty(H_hr.HcoeL)
            H_hr.HcoeL = conj(H_hr.HcoeL(:,:,H_hr.Duality_vector_dist));
        end
        if ~isempty(H_hr.HnumL)
            H_hr.HnumL = conj(H_hr.HnumL(:,:,H_hr.Duality_vector_dist));
        end
        
    elseif strcmp(H_hr.Type, 'sparse')
        % Sparse matrix conjugation
        if ~isempty(H_hr.HcoeL)
            H_hr.HcoeL = cellfun(@conj, H_hr.HcoeL(H_hr.Duality_vector_dist), 'UniformOutput', false);
        end
        if ~isempty(H_hr.HnumL)
            H_hr.HnumL = cellfun(@conj, H_hr.HnumL(H_hr.Duality_vector_dist), 'UniformOutput', false);
        end
        
    else
        error('Unsupported matrix type: %s', H_hr.Type);
    end
end