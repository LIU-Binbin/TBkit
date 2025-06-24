classdef pauli_matrix
    % Pauli matrix class representing the 4 Pauli matrices
    % and supporting common matrix operations like addition,
    % multiplication, inversion, etc.
    
    properties
        mat;  % Matrix representing a Pauli matrix
    end
    
    methods
        % Constructor for Pauli matrix
        function PM = pauli_matrix(label)
            if nargin < 1
                % Default constructor: creates an array of Pauli matrices (I, X, Y, Z)
                for i = 1:4
                    PM(i) = pauli_matrix(i - 1);
                end
            else
                if isa(label,'pauli_matrix')
                    PM.mat = label.mat;
                    return;
                end
                % Initialize matrix based on the given label
                if isscalar(label)
                    switch label
                        case {'0', 'I', 0}
                            PM.mat = [1, 0; 0, 1];  % Identity matrix
                        case {'x', '1', 'X', 1}
                            PM.mat = [0, 1; 1, 0];  % Pauli-X matrix
                        case {'y', '2', 'Y', 2}
                            PM.mat = [0, -1i; 1i, 0];  % Pauli-Y matrix
                        case {'z', '3', 'Z', 3}
                            PM.mat = [1, 0; 0, -1];  % Pauli-Z matrix
                        otherwise

                            error('Invalid label for Pauli matrix');
                    end
                else
                    PM.mat = label;
                end
            end
        end
    end
    
    %% Operator overloads for common matrix operations
    methods
        % Display method for printing matrix
        function disp(PM)
            for i = 1:length(PM)
                disp(PM(i).mat);
            end
        end
        
        % Overload addition operator
        function PM = plus(PM1, PM2)
            if isa(PM1, 'pauli_matrix') && isa(PM2, 'pauli_matrix')
                PM = PM1;
                PM.mat = PM1.mat + PM2.mat;
            else
                error('Invalid operands for + operator');
            end
        end
        
        % Overload subtraction operator
        function PM = minus(PM1, PM2)
            if isa(PM1, 'pauli_matrix') && isa(PM2, 'pauli_matrix')
                PM = PM1;
                PM.mat = PM1.mat - PM2.mat;
            else
                error('Invalid operands for - operator');
            end
        end
        
        % Unary minus operator (negation)
        function PM = uminus(PM)
            PM.mat = -PM.mat;
        end
        
        % Overload element-wise multiplication operator
        function PM = times(PM1, PM2)
            if isa(PM1, 'pauli_matrix') && isa(PM2, 'pauli_matrix')
                PM = PM1;
                PM.mat = PM1.mat * PM2.mat;
            else
                error('Invalid operands for .* operator');
            end
        end
        
        % Overload matrix multiplication operator
        function PM = mtimes(PM1, PM2)
            if isa(PM1, 'pauli_matrix') && isa(PM2, 'pauli_matrix')
                PM = PM1;
                PM.mat = kron(PM1.mat, PM2.mat);  % Kronecker product
            elseif ~isa(PM1, 'pauli_matrix') && isa(PM2, 'pauli_matrix')
                PM = PM2;
                PM.mat = (PM1 * PM2.mat);
            elseif   isa(PM1, 'pauli_matrix') && ~isa(PM2, 'pauli_matrix')
                PM = PM1;
                PM.mat =(PM1.mat * PM2);
            else
                error('Invalid operands for * operator');
            end
        end
        
        % Overload right matrix division operator
        function PM = mrdivide(PM1, PM2)
            if isa(PM1, 'pauli_matrix') && isa(PM2, 'pauli_matrix')
                PM = PM1;
                PM.mat = PM1.mat / PM2.mat;
            else
                error('Invalid operands for / operator');
            end
        end
        
        % Convert to double precision matrix
        function mat = double(PM)
            mat = PM.mat;
        end
        
        % Return symbolic representation of the matrix
        function symmat = sym(PM)
            if length(PM) == 1
                symmat = PM.mat;
            else
                for i=1:length(PM)
                    symmat(:,:,i) = PM(i).mat;
                end
            end
            symmat = sym(symmat);
        end
        % Compute inverse of the matrix
        function PM = inv(PM)
            PM.mat = inv(PM.mat);
        end
        
        % Compute matrix exponential
        function PM = expm(PM)
            PM.mat = expm(PM.mat);
        end
        
        % Compute Hermitian transpose (conjugate transpose)
        function PM = ctranspose(PM)
            PM.mat = PM.mat';
        end
        
        % Compute transpose of the matrix
        function PM = transpose(PM)
            PM.mat = PM.mat.';
        end
        
        % Compute the anticommutator [PM1, PM2] = PM1*PM2 + PM2*PM1
        function PM = anticommute(PM1, PM2)
            PM = PM1;
            PM.mat = PM1.mat * PM2.mat + PM2.mat * PM1.mat;
        end
        
        % Compute the commutator [PM1, PM2] = PM1*PM2 - PM2*PM1
        function PM = commute(PM1, PM2)
            PM = PM1;
            PM.mat = PM1.mat * PM2.mat - PM2.mat * PM1.mat;
        end
    end
    
    %% Static methods for creating symbolic representations
    methods(Static)
        % Generate matrix for S transformation (related to Pauli matrices)
        function Smat_inv = S()
            Smat = zeros(4);
            Pauli_L = sym(pauli_matrix);
            for i = 1:4
                tmp_mat = Pauli_L(:,:,i);
                tmp_mat_r = real(tmp_mat);
                tmp_mat_i = imag(tmp_mat);
                Smat(i,1:2) = diag(tmp_mat_r);
                Smat(i,3) = diag(tmp_mat_r, 1);
                Smat(i,4) = diag(tmp_mat_i, 1);
            end
            Smat_inv = inv(Smat);
        end
        
        % Generate symbolic list of Pauli matrices
        function Pauli_L = L()
            count = 0;
            Pauli_L = sym(zeros(1, 4));
            for i = 0:3
                count = count + 1;
                Pauli_L(count) = sym(['Pauli_', num2str(i)]);
            end
        end
    end
end
