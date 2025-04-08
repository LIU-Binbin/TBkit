classdef Term
    % The Term class is a MATLAB class designed to represent a mathematical term that consists of two key components:
    % 
    % Symbolic Polynomial (symbolic_polynomial): A symbolic expression that represents the mathematical part of the term, typically involving variables and constants. The symbolic polynomial can be manipulated symbolically for further analysis.
    % 
    % Pauli Matrix (pauli_mat): A 4x4(n*n) matrix associated with the term, which is typically a Pauli matrix or related to some quantum mechanical operator. Pauli matrices are commonly used in quantum mechanics and quantum computing to represent spin operators or other observables in a two-level quantum system.
    properties
        symbolic_polynomial; % Symbolic polynomial expression for the term
        pauli_mat;          % Pauli matrix associated with the term (4x4 matrix)
    end
    
    methods
        % Constructor to initialize a Term object
        function term = Term(symbolic_polynomial, pauli_mat)
            if nargin < 1
                % If no input arguments, initialize with default values
                term.symbolic_polynomial = sym(0); % Default symbolic polynomial (0)
                term.pauli_mat = zeros(4);         % Default Pauli matrix (4x4 zero matrix)
            else
                % Otherwise, initialize with provided values
                term.symbolic_polynomial = symbolic_polynomial;
                term.pauli_mat = pauli_mat;
            end
        end

        % Display method for Term object
        function disp(term)
            % Display symbolic polynomial and Pauli matrix for the term
            disp('Symbolic Polynomial:');
            disp(term.symbolic_polynomial);
            disp('Pauli Matrix:');
            disp(term.pauli_mat);
        end

        % Method to generate the symbolic matrix representation of the term
        function symmat = sym(term)
            % Initialize a symbolic matrix to accumulate contributions
            symmat = sym(zeros(size(term(1).pauli_mat))); % Assuming Pauli matrix is square
            for i = 1:length(term)
                % Sum each Term's symbolic polynomial multiplied by its Pauli matrix
                symmat = symmat + term(i).symbolic_polynomial * term(i).pauli_mat;
            end
        end
        
        % Overload the plus operator to add two Term objects
        function C = plus(A, B)
            if isa(A, 'Term') && isa(B, 'Term')
                % If both A and B are Term objects, concatenate them
                C = [A, B];
            elseif isa(B, 'Term') && ~isa(A, 'Term')
                % If B is a Term and A is not, handle accordingly (e.g., A might be HK)
                C = A + B;
            elseif ~isa(B, 'Term') && isa(A, 'Term')
                % If A is a Term and B is not, handle accordingly
                C = A + B;
            else
                % Return zero if neither A nor B is a Term object
                C = 0;
            end
        end

        % Compute the commutator (or anti-commutator) of two Pauli matrices
        function result = commute(A, B, isanti)
            % Set default value for isanti if not provided (0 for commutator, 1 for anti-commutator)
            if nargin < 3
                isanti = 0; % Default is commutator
            end
            
            % If A or B is a Pauli matrix object, extract its matrix
            if isa(A, 'pauli_matrix') 
                A = A.mat;
            end
            if isa(B, 'pauli_matrix')
                B = B.mat;
            end
            
            % Convert matrices to symbolic for easier manipulation
            M1 = sym(A);
            M2 = sym(B);
            
            % Compute the commutator or anti-commutator based on the flag `isanti`
            if isanti == 1
                % Anti-commutator: (B*A) - (A*B)
                msymtmp = (M2 * conj(M1)) - (M1 * M2);
            else
                % Commutator: (B*A) - (A*B)
                msymtmp = (M2 * M1) - (M1 * M2);
            end
            
            % Normalize the resulting matrix by dividing by the first element
            mattmp = msymtmp / msymtmp(1,1);
            tf = diag(mattmp); % Check the diagonal elements
            
            % Check if all diagonal elements are 1 (identity-like)
            label = all(tf == ones(length(tf), 1));
            
            % If the result is an identity-like matrix, return its eigenvalue
            if label
                eigentmp = diag(msymtmp);
                result = eigentmp(1);
            else
                % If not, return the commutator or anti-commutator expressions
                disp('Dependent')
                if isanti == 1
                    % Anti-commutator: (B*A) - (A*B) and (B*A) + (A*B)
                    result{1} = (M2 * conj(M1)) - (M1 * M2);
                    result{2} = (M2 * conj(M1)) + (M1 * M2);
                else
                    % Commutator: (B*A) - (A*B) and (B*A) + (A*B)
                    result{1} = (M2 * M1) - (M1 * M2);
                    result{2} = (M2 * M1) + (M1 * M2);
                end
            end
        end

        % Method to shift terms based on a given k-point (translation in k-space)
        function term_list = subsk(term_list, kpoints_r)
            % Define symbolic variables for k-space
            syms k_x k_y k_z;
            
            % Shift the k-points by the provided values
            k_x = k_x - kpoints_r(1);
            k_y = k_y - kpoints_r(2);
            k_z = k_z - kpoints_r(3);
            
            % Substitute the new k-values into each term's symbolic polynomial
            for i = 1:length(term_list)
                term_list(i).symbolic_polynomial = subs(term_list(i).symbolic_polynomial, ...
                                                        {k_x, k_y, k_z}, ...
                                                        {kpoints_r(1), kpoints_r(2), kpoints_r(3)});
            end
        end
    end
end
