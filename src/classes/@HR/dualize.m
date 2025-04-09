function H_hr = dualize(H_hr)
%DUALIZE Generate dual Hamiltonian by reversing lattice vectors and establishing duality mapping
% This function processes the input H_hr object to create a dual structure by:
% 1. Generating opposite vectors for all lattice vectors
% 2. Handling special column swapping for 5-column format vectors
% 3. Creating duality mapping between original and dual vectors
% 4. Expanding vector list with missing dual vectors when necessary
%
% Input: H_hr - HR Hamiltonian object containing:
% NRPTS: Number of lattice vectors
% vectorL: List of lattice vectors [NxDim or Nx5 matrix]
% Dim: Dimension of vectors
% Output: Modified H_hr with added dual vectors and Duality_vector_dist property

NRPTS_ = H_hr.NRPTS; % Get original number of lattice vectors
vectorList = H_hr.vectorL; % Extract original vector list

% Create opposite vectors by negating spatial components
vectorList_oppo(:,1:H_hr.Dim) = -vectorList(:,1:H_hr.Dim);

% Handle special 5-column format (3 spatial + 2 extra dimensions)
if size(vectorList,2) == 5
% Swap last two columns for opposite vectors
vectorList_oppo(:,H_hr.Dim+1) = vectorList(:,H_hr.Dim+2);
vectorList_oppo(:,H_hr.Dim+2) = vectorList(:,H_hr.Dim+1);
end

% Initialize duality mapping vector
H_hr.Duality_vector_dist = zeros(NRPTS_,1);

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
    H_hr.Duality_vector_dist(j) = i;            % Record reverse mapping
end

% Store duality mapping: original -> dual
H_hr.Duality_vector_dist(i) = j;
end
end