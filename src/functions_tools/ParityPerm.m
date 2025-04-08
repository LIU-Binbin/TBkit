function Parity = ParityPerm(permVec)
    % ParityPerm Calculates the parity (sign) of a permutation
    %   Input: permVec - a vector representing the permutation
    %   Output: Parity - +1 for even permutations, -1 for odd permutations

    % Check if input is a matrix and convert to a vector if necessary
    if ~isvector(permVec)
        permVec = permMatToVec(permVec);  % Convert permutation matrix to vector if input is a matrix
    end
    
    % Initialize parity as +1 (even permutation)
    Parity = 1;

    % Count inversions efficiently
    n = numel(permVec);
    for i = 1:n-1
        for j = i+1:n
            if permVec(i) > permVec(j)
                Parity = -Parity;  % Flip parity for each inversion
            end
        end
    end
end
