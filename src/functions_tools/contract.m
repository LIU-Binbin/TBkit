function [T] = contract(A, I1, B, I2, varargin)
% CONTRACT performs the contraction of two tensors A and B over the specified indices I1 and I2.
% It supports partial or full contraction, and the resulting tensor can be permuted.
%
% Input:
%   A       - Tensor A.
%   I1      - Indices for contraction on tensor A (can be an array).
%   B       - Tensor B.
%   I2      - Indices for contraction on tensor B (must match the size of I1).
%   varargin - Optional argument specifying the desired permutation of the resulting tensor.
%
% Output:
%   T       - Resulting contracted tensor (with the desired permutation if provided).

    % Handle the case when both I1 and I2 are empty, implying full contraction.
    if isempty(I1) && isempty(I2)
        I1 = 1:numel(size(A));
        I2 = 1:numel(size(B));
    end

    % Ensure that the contraction indices are valid (i.e., they must match in size).
    if numel(I1) ~= numel(I2)
        error('The lengths of I1 and I2 must be equal!');
    end

    % Get the dimensions of tensors A and B.
    DimA = size(A);
    DimB = size(B);

    % Check if the dimensions for contraction match.
    if isequal(DimA(I1), DimB(I2))
        % Create sorted index arrays for A and B after contraction.
        Ia = setdiff(1:numel(DimA), I1);
        Ib = setdiff(1:numel(DimB), I2);

        % Ensure that I1 indices are part of Ia, and I2 indices are part of Ib.
        if any(~ismember(I1, Ia)) || any(~ismember(I2, Ib))
            warning('Some indices for contraction are not valid!');
        end

        % Permute A and B according to the specified contraction indices.
        A = permute(A, [Ia, I1]);
        B = permute(B, [I2, Ib]);

        % Reshape the tensors to perform matrix multiplication.
        A = reshape(A, [], prod(DimA(I1)));
        B = reshape(B, prod(DimB(I2)), []);

        % Perform the contraction as matrix multiplication.
        t = A * B;

        % Get the dimensions of the resulting tensor after contraction.
        DimAR = DimA(Ia);
        DimBR = DimB(Ib);
        
        % Reshape the result back to a tensor if necessary.
        if ~isempty([DimAR, DimBR])
            T = reshape(t, [DimAR, DimBR]);
        else
            T = t;  % If the result is just a scalar, return the scalar.
        end

        % Apply any optional permutation to the resulting tensor.
        if ~isempty(varargin)
            T = permute(T, varargin{1});
        end

    else
        error('Contraction indices do not match between tensors A and B!');
    end
end

