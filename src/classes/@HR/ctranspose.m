function H_hr = ctranspose(H_hr)
% CTRANSPOSE Hermitian conjugate transpose of HR object
%
%   H_conj = CTRANSPOSE(H_hr) computes the conjugate transpose of the 
%   Hamiltonian that satisfies the Hermitian symmetry condition:
%       H_conj(r) = H_hr†(-r)
%
%   Features:
%   - Handles both vector hopping and matrix representations
%   - Supports list and matrix storage formats
%   - Maintains internal consistency of HR object properties
%
%   Input:
%       H_hr - HR object with required properties:
%           dualize() method, Duality_vector_dist permutation vector
%           vectorhopping flag, AvectorL/BvectorL/CvectorL (if vectorhopping)
%           Type (storage format), num/coe flags, HnumL/HcoeL
%
%   Output:
%       H_hr - Conjugate transposed HR object
%
%   Example:
%       H_conj = H_hr';  % Equivalent to ctranspose(H_hr)
%
%   See also DUALIZE, CONJ, AUTOHERMI

    % Validate dualization prerequisites
    try
        H_hr = H_hr.dualize();
    catch ME
        error('CTRANSPOSE failed during dualization: %s', ME.message);
    end
    
    % Check required permutation vector
    if ~isprop(H_hr, 'Duality_vector_dist') || isempty(H_hr.Duality_vector_dist)
        error('Missing Duality_vector_dist property');
    end

    % Main conjugation operations
    if H_hr.vectorhopping
        process_vector_hopping();
    else
        process_matrix_hopping();
    end

    % Nested processing functions
    function process_vector_hopping()
        % Validate vector hopping components
        validate_vector_components();
        
        % Process A/B vectors
        H_hr.AvectorL = H_hr.AvectorL(H_hr.Duality_vector_dist, :);
        H_hr.BvectorL = -H_hr.BvectorL(H_hr.Duality_vector_dist, :);
        
        % Process C vectors with complex conjugation
        process_c_vectors();
    end

    function process_matrix_hopping()
        switch H_hr.Type
            case 'list'
                process_list_format();
            otherwise
                process_matrix_format();
        end
    end

    % Helper functions
    function validate_vector_components()
        required_fields = {'AvectorL', 'BvectorL', 'CvectorL'};
        for f = required_fields
            if ~isfield(H_hr, f{1}) || isempty(H_hr.(f{1}))
                error('Missing vector component: %s', f{1});
            end
        end
        
        % Validate CvectorL structure
        if mod(size(H_hr.CvectorL, 1), 2) ~= 0
            error('CvectorL must have even number of rows for real/imag separation');
        end
    end

    function process_c_vectors()
        % Split and process complex components
        mid_idx = size(H_hr.CvectorL, 1)/2;
        CLr = H_hr.CvectorL(1:mid_idx, :);
        CLi = H_hr.CvectorL(mid_idx+1:end, :);
        
        % Apply permutation and conjugation
        H_hr.CvectorL = [
            CLr(H_hr.Duality_vector_dist, :); 
            -CLi(H_hr.Duality_vector_dist, :)
        ];
    end

    function process_list_format()
        % Process list-format Hamiltonian
        if H_hr.num
            H_hr.HnumL = conj(H_hr.HnumL(H_hr.Duality_vector_dist));
        end
        if H_hr.coe
            H_hr.HcoeL = conj(H_hr.HcoeL(H_hr.Duality_vector_dist));
        end
    end

    function process_matrix_format()
        % Process matrix-format Hamiltonian
         % NRPTS_ = H_hr.NRPTS; % Get original number of lattice vectors
        if H_hr.coe
            tmpHcoeL = H_hr.HcoeL;
            for i = 1:length(H_hr.Duality_vector_dist)
                tmpHcoeL(:,:,i)= H_hr.HcoeL(:, :, H_hr.Duality_vector_dist(i))';
            end
            H_hr.HcoeL = tmpHcoeL;
        end
        if H_hr.num
            tmpHnumL = H_hr.HnumL;
            for i = 1:length(H_hr.Duality_vector_dist)
                tmpHnumL(:,:,i)= H_hr.HnumL(:, :, H_hr.Duality_vector_dist(i))';
            end
            H_hr.HnumL = tmpHnumL;
        end
    end
end