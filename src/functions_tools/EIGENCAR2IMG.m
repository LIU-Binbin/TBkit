function IMG = EIGENCAR2IMG(EIGENCAR, dE, Erange, kselect)
    % EIGENCAR2IMG - Convert eigenvalue data to an image representation
    %
    % Syntax: IMG = EIGENCAR2IMG(EIGENCAR, dE, Erange, kselect)
    %
    % Inputs:
    %   EIGENCAR   - Matrix of eigenvalues (Bands x k-points)
    %   dE         - Energy step size (default: 0.1)
    %   Erange     - Energy range [Emin, Emax] (default: min and max of EIGENCAR)
    %   kselect    - Indices of k-points to select (default: all columns)
    %
    % Outputs:
    %   IMG        - 2D image representation of eigenvalue distribution
    %
    % Example:
    %   IMG = EIGENCAR2IMG(EIGENCAR, 0.1, [0,10], 1:5);
    %

    % Set default values for optional parameters
    if nargin < 2
        dE = 0.1;
    end
    if nargin < 3
        Erange = [min(EIGENCAR(:)) max(EIGENCAR(:))];
    end
    if nargin < 4
        kselect = 1:size(EIGENCAR, 2);
    end

    % Get dimensions of the eigenvalue data
    [NBANDS, NKPTS] = size(EIGENCAR);
    DATA = EIGENCAR(:, kselect);
    X_nodes = size(DATA, 2);
    Y_nodes = round(abs(Erange(2) - Erange(1)) / dE);
    Emin = Erange(1);

    % Initialize sparse matrix for efficient data accumulation
    IMG_sparse = sparse(Y_nodes, X_nodes);

    % Populate the sparse matrix
    for i = 1:X_nodes
        for j = 1:NBANDS
            energy = DATA(j, i);
            NE = round((energy - Emin) / dE);
            if NE >= 1 && NE <= Y_nodes
                IMG_sparse(NE, i) = IMG_sparse(NE, i) + 0.1;
            end
        end
    end

    % Convert sparse matrix to full matrix for output
    IMG = full(IMG_sparse);
end



