function KronVector = SemiProductVector(VectorLeft, VectorRight)
    % Compute the Kronecker product of VectorLeft with a column vector of ones
    % to expand VectorLeft.
    expandedLeft = kron(VectorLeft, ones(size(VectorRight, 1), 1));
    
    % Repeat VectorRight to match the expanded size of VectorLeft.
    repeatedRight = repmat(VectorRight, [size(VectorLeft, 1), 1]);
    
    % Concatenate the expanded VectorLeft and the repeated VectorRight.
    KronVector = [expandedLeft, repeatedRight];
end
