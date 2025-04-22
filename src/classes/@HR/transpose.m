function H_hr = transpose(H_hr)
%TRANSPOSE Overloaded transpose for HR objects
%
%   H_HR = TRANSPOSE(H_HR) computes the non-conjugate transpose
%   of the Hamiltonian matrices in an HR object.
%
%   Input:
%       H_hr - HR object to transpose
%
%   Output:
%       H_hr - Transposed HR object
%
%   See also HR, CTRANSPOSE, PAGETRANSPOSE
if strcmp(H_hr.Type,'mat')
    if H_hr.coe
        H_hr.HcoeL = pagetranspose(H_hr.HcoeL);
    end
    if H_hr.num
        H_hr.HnumL = pagetranspose(H_hr.HnumL);
    end
elseif strcmp(H_hr.Type,'list')
    NRPTS_ = H_hr.NRPTS; % Get original number of lattice vectors
    vectorList = H_hr.vectorL; % Extract original vector list
    % Create opposite vectors by negating spatial components
    vectorList_oppo(:,1:H_hr.Dim) = vectorList(:,1:H_hr.Dim);

    % Handle special 5-column format (3 spatial + 2 extra dimensions)
    if size(vectorList,2) == 5
        % Swap last two columns for opposite vectors
        vectorList_oppo(:,H_hr.Dim+1) = vectorList(:,H_hr.Dim+2);
        vectorList_oppo(:,H_hr.Dim+2) = vectorList(:,H_hr.Dim+1);
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
