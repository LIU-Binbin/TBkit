function [Asort, Usort] = sorteig(U, A)
    % Determine the mode based on the number of input arguments
    if nargin < 2
        mode = 'eigenval';
        NUM_WAN = length(U);
    else
        mode = 'whole';
        NUM_WAN = length(A);
    end
    NBANDS = length(U);
    
    % Check if U is a vector or a matrix
    if ~isvector(U)
        SortTmp = diag(U); % Extract eigenvalues from the diagonal
        vec = false;
    else
        SortTmp = U;
        vec = true;
    end
    
    % Sorting based on the mode
    if strcmp(mode, 'whole')
        if size(U, 2) ~= size(A, 2)
            % Handle non-Hermitian matrices
            [Usort, IJ] = sort(SortTmp, 1, 'ComparisonMethod', 'real');
            if ~vec
                Usort = diag(Usort);
            end
            Asort = zeros(NUM_WAN, NBANDS);
        else
            % Sort eigenvalues and rearrange eigenvectors
            Asort = zeros(NUM_WAN, NBANDS);
            [Usort, IJ] = sort(SortTmp, 1, 'ComparisonMethod', 'real');
            for jj = 1:NBANDS
                Asort(:, jj) = A(:, IJ(jj));
            end
            if ~vec
                Usort = diag(Usort);
            end
        end
    elseif strcmp(mode, 'eigenval')
        % Sort eigenvalues in descending order
        SortTmp = diag(U);
        [Usort, ~] = sort(SortTmp, 1, 'ComparisonMethod', 'real');
        if ~vec
            Usort = diag(Usort);
        end
        Asort = [];
    end
end
