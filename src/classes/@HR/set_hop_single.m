function H_hr = set_hop_single(H_hr, amp, hi, hj, vector, mode)
%SET_HOP_SINGLE Set hopping parameters in tight-binding Hamiltonian
%   This function modifies the hopping terms in a tight-binding Hamiltonian
%   structure (H_hr) based on specified parameters and operation mode.
%
%   Inputs:
%   H_hr   - Hamiltonian structure containing:
%            Type: Data storage type ('list'/'sparse'/'matrix')
%            vectorL: List of R vectors
%            NRPTS: Number of R points
%            HnumL: Numerical hopping parameters
%            HcoeL: Symbolic hopping parameters
%   amp    - Hopping amplitude (numeric or symbolic)
%   hi,hj  - Orbital indices for hopping term
%   vector - Lattice vector for hopping term [n1,n2,...]
%   mode   - Operation mode:
%            'set'    : Overwrite existing value
%            'add'    : Add to existing value
%            'sym'    : Set symbolic value
%            'symadd' : Add symbolic value
%
%   Output:
%   H_hr   - Updated Hamiltonian structure
%
%   Features:
%   - Handles both numeric and symbolic coefficients
%   - Supports multiple storage formats (list/sparse/matrix)
%   - Automatic expansion of R-point list when new vectors are added

% Process input vector based on storage type
V = H_hr.vectorL;
if strcmp(H_hr.Type, 'list')
    vector = [vector, hi, hj];  % Combine indices with lattice vector
end

% Find or create R-point entry
if isempty(V)
    % Initialize first R-point entry
    seq = H_hr.NRPTS + 1;
    H_hr = H_hr.add_empty_one(vector);
else
    % Check existing R-points
    [~, seq] = ismember(vector, V, 'rows');

    % Create new entry if vector not found
    if seq == 0
        seq = H_hr.NRPTS + 1;
        H_hr = H_hr.add_empty_one(vector);
    end
end

% Update hopping parameters based on storage type and mode
switch H_hr.Type
    case 'list'
        handle_list_storage();
    case 'sparse'
        handle_sparse_storage();
    otherwise  % Matrix storage
        handle_matrix_storage();
end

% Nested functions for different storage types
    function handle_list_storage()
        switch mode
            case 'set'
                check_overwrite_warning(@() H_hr.HnumL(seq), amp);
                H_hr.HnumL(seq) = amp;
            case 'add'
                H_hr.HnumL(seq) = H_hr.HnumL(seq) + amp;
            case 'sym'
                check_overwrite_warning(@() H_hr.HcoeL(seq), amp);
                H_hr.HcoeL(seq) = amp;
            case 'symadd'
                H_hr.HcoeL(seq) = H_hr.HcoeL(seq) + amp;
        end
    end

    function handle_sparse_storage()
        switch mode
            case 'set'
                check_overwrite_warning(@() H_hr.HnumL{seq}(hi,hj), amp);
                H_hr.HnumL{seq}(hi,hj) = amp;
            case 'add'
                H_hr.HnumL{seq}(hi,hj) = H_hr.HnumL{seq}(hi,hj) + amp;
            case 'sym'
                check_overwrite_warning(@() H_hr.HcoeL{seq}(hi,hj), amp);
                H_hr.HcoeL{seq}(hi,hj) = amp;
            case 'symadd'
                H_hr.HcoeL{seq}(hi,hj) = H_hr.HcoeL{seq}(hi,hj) + amp;
        end
    end

    function handle_matrix_storage()
        switch mode
            case 'set'
                check_overwrite_warning(@() H_hr.HnumL(hi,hj,seq), amp);
                H_hr.HnumL(hi,hj,seq) = amp;
            case 'add'
                H_hr.HnumL(hi,hj,seq) = H_hr.HnumL(hi,hj,seq) + amp;
            case 'sym'
                check_overwrite_warning(@() H_hr.HcoeL(hi,hj,seq), amp);
                H_hr.HcoeL(hi,hj,seq) = amp;
            case 'symadd'
                H_hr.HcoeL(hi,hj,seq) = H_hr.HcoeL(hi,hj,seq) + amp;
        end
    end

% Helper function for overwrite warnings
    function check_overwrite_warning(get_value_fn, new_value)
        if get_value_fn() ~= 0 && ~isequal(get_value_fn(), new_value)
            warning('Potential overwrite detected. Consider using add mode.');
            fprintf('Conflict at vector: %s\n', num2str(vector));
        end
    end

end