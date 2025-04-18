function mat = pauli_matrice(dir)
    % pauli_matrixe - Generates Pauli matrices or specific matrix based on input
    % 
    % Usage:
    %   sigma = pauli_matrixe();            % Returns all Pauli matrices as a struct
    %   sigma_i = pauli_matrixe('x');       % Returns specific Pauli matrix based on input
    %
    % Input:
    %   dir (optional) - Specifies which Pauli matrix to return.
    %                    Options: '0' or 'o' for identity matrix (sigma_0)
    %                            '1' or 'x' for Pauli X matrix (sigma_x)
    %                            '2' or 'y' for Pauli Y matrix (sigma_y)
    %                            '3' or 'z' for Pauli Z matrix (sigma_z)
    %
    % Output:
    %   mat - The requested Pauli matrix or struct containing all matrices.
    %
    % Example:
    %   tau = pauli_matrixe();
    %   tau.o;  % Identity matrix
    %   tau.x;  % Pauli X matrix
    %   tau.y;  % Pauli Y matrix
    %   tau.z;  % Pauli Z matrix
    %
    % Notes:
    %   - If no input is provided, the function returns all four matrices as a struct.
    %   - Valid inputs are '0', 'o', '1', 'x', '2', 'y', '3', and 'z'.
    
    % Imaginary unit
    I = 1i;
    
    % If no argument is provided, return all Pauli matrices in a struct
    if nargin < 1
        sigma.o = eye(2);                         % Identity matrix
        sigma.x = [0 1; 1 0];                     % Pauli X
        sigma.y = [0 -I; I 0];                    % Pauli Y
        sigma.z = [1 0; 0 -1];                    % Pauli Z
        mat = sigma;
        return;
    end
    
    % Check the input and return the corresponding Pauli matrix
    switch lower(char(string(dir)))
        case {'0', 'o'}
            mat = eye(2);                           % Identity matrix
        case {'1', 'x'}
            mat = [0 1; 1 0];                      % Pauli X
        case {'2', 'y'}
            mat = [0 -I; I 0];                     % Pauli Y
        case {'3', 'z'}
            mat = [1 0; 0 -1];                     % Pauli Z
        otherwise
            error('Invalid input. Valid inputs are ''0'', ''1'', ''2'', ''3'', or ''o'', ''x'', ''y'', ''z''.');
    end
end
