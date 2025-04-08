classdef Trig
    % The Trig class is a MATLAB class designed to represent a mathematical term that consists of two key components:
    % 
    % Symbolic Polynomial (symbolic_polynomial): trigonometric function
    % 
    % Pauli Matrix (pauli_mat): A 4x4(n*n) matrix associated with the term, which is typically a Pauli matrix or related to some quantum mechanical operator. Pauli matrices are commonly used in quantum mechanics and quantum computing to represent spin operators or other observables in a two-level quantum system.
    
    properties
        symbolic_polynomial;  % Symbolic polynomial representation (e.g., trigonometric function)
        pauli_mat;            % Pauli matrix or other matrix representation
    end
    
    methods
        % Constructor: Initialize the symbolic polynomial and Pauli matrix
        function Trig = Trig(symbolic_polynomial, pauli_mat)
            % Default constructor initializes to zero polynomial and NaN for matrix if no arguments
            if nargin < 1
                Trig.symbolic_polynomial = sym(0);  % Default symbolic polynomial (zero)
                Trig.pauli_mat = NaN;               % Default matrix (NaN indicating undefined)
                Trig(1) = [];                       % Empty array if no input arguments
            else
                % Assign provided symbolic polynomial and Pauli matrix
                Trig.symbolic_polynomial = symbolic_polynomial;
                Trig.pauli_mat = pauli_mat;
            end
        end
        
        % Unary minus operator: Negates the symbolic polynomial
        function A = uminus(A)
            % Negates the symbolic polynomial of the object
            A.symbolic_polynomial = -A.symbolic_polynomial;
        end
        
        % Addition operator: Adds two Trig objects or a Trig object and another type (e.g., HK)
        function C = plus(A, B)
            % Case 1: Adding two Trig objects
            if isa(A, 'Trig') && isa(B, 'Trig')
                C = A;  % Copy A to C
                % Append elements from B to the result
                C(length(A)+1 : length(A)+length(B)) = B;
                
            % Case 2: Adding a Trig object to something else (e.g., HK)
            elseif isa(B, 'Trig') && ~isa(A, 'Trig')
                disp('Invalid addition: Handling case with non-Trig type A');
                if isa(A, 'HK')  % Special case for HK (if applicable)
                    C = A.plus(B);
                else
                    C = A;  % No addition, just return A as is
                end
                
            % Case 3: Adding something else (e.g., HK) to a Trig object
            elseif ~isa(B, 'Trig') && isa(A, 'Trig')
                if isa(B, 'HK')  % Special case for HK (if applicable)
                    if length(B) == 2
                        C = B.plus(A);
                    else
                        C = B;  % Just return B if not exactly 2
                    end
                else
                    C = B;  % Return B as is
                end
            else
                % Case where neither A nor B are of the expected types
                C = 0;  % Return zero as a fallback
            end
        end
    end
end
