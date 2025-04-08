function permVec = permMatToVec(permMat)
    % permMatToVec Converts a permutation matrix to a permutation vector
    %   Input: permMat - a square permutation matrix
    %   Output: permVec - a vector representing the permutation
    % Define a permutation matrix
    % Example:
    %   permMat = [
    %     0 1 0;
    %     0 0 1;
    %     1 0 0
    %     ];
    % 
    % % Call the permMatToVec function to convert the matrix to a permutation vector
    %   permVec = permMatToVec(permMat);
    % 
    % % Display the result
    %   disp('Permutation vector:');
    %   disp(permVec);
    [n, m] = size(permMat);  % Get the size of the matrix (assuming square matrix)
    
    if n ~= m
        error('Input matrix must be square');
    end

    % Find the column index (position) of the 1 in each row (this gives the permutation vector)
    permVec = find(permMat.' == 1);  % Use transpose to find column indices
end        