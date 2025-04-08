function permMat = permVecToMat(permVec)
    % permVecToMat Converts a permutation vector to a permutation matrix
    %   Input: permVec - a vector of integers representing the permutation
    %   Output: permMat - a square matrix with ones placed according to the permutation

    n = length(permVec);  % Get the length of the permutation vector
    permMat = eye(n);     % Initialize the matrix as an identity matrix

    % Update the identity matrix to match the permutation
    permMat = permMat(:, permVec);
end

