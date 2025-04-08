classdef quaternion
    % This class represents a quaternion and provides various operations
    % such as addition, multiplication, conjugation, and inversion.
    
    properties
        U  % The quaternion components (real and imaginary parts)
    end
    
    properties(Dependent)
        expression  % Dependent property to return the quaternion as a symbolic expression
    end
    
    methods
        % Constructor: Initialize quaternion from real and imaginary parts
        function Q = quaternion(Re, I, J, K)
            % If no arguments are provided, set a default quaternion
            if nargin < 1
                Q.U = [0; 0; 0; 0];  % Default quaternion [0; 0; 0; 0]
            else
                Q.U = [Re; I; J; K];  % Initialize with provided components
            end
        end
    end
    
    % Getter method for 'expression' dependent property
    methods
        function expression = get.expression(Q)
            % Returns the quaternion in symbolic form as a polynomial
            syms I J K;
            expression = Q.U(1) + Q.U(2)*I + Q.U(3)*J + Q.U(4)*K;
        end
    end
    
    % Static methods to generate standard quaternion units
    methods(Static)
        function Q_I = I()
            % Returns the quaternion representing the imaginary unit i
            Q_I = quaternion(0, 1, 0, 0);
        end
        
        function Q_J = J()
            % Returns the quaternion representing the imaginary unit j
            Q_J = quaternion(0, 0, 1, 0);
        end
        
        function Q_K = K()
            % Returns the quaternion representing the imaginary unit k
            Q_K = quaternion(0, 0, 0, 1);
        end
    end
    
    % Methods for quaternion operations
    methods
        % Conjugate: Negate the imaginary components (i, j, k)
        function Q = conj(Q)
            Q.U(2:4) = -Q.U(2:4);  % Conjugate by negating the imaginary parts
        end
        
        % Inner product: Compute the norm of the product of two quaternions
        function value = innerproduct(Q1, Q2)
            value = norm(Q1 * conj(Q2));  % Return the norm of the product
        end
        
        % Norm: Return the magnitude of the quaternion
        function normValue = norm(Q)
            normValue = norm(Q.U);  % Compute the norm of the quaternion
        end
        
        % Inverse: Return the inverse of the quaternion
        function Q_inv = inverse(Q)
            Q_inv = 1 / norm(Q)^2 * conj(Q);  % Inverse is conjugate divided by the norm squared
        end
        
        % Real part: Return the real component of the quaternion
        function Re = real(Q)
            Re = Q.U(1);  % The real part is the first component
        end
        
        % Quaternion multiplication (overloaded '*')
        function Q = mtimes(Q1, Q2)
            % Perform quaternion multiplication
            if isa(Q1, 'quaternion') && isa(Q2, 'quaternion')
                % Quaternion * Quaternion multiplication
                q1 = Q1.U;
                q2 = Q2.U;
                
                Q.U(1) = q1(1) * q2(1) - dot(q1(2:4), q2(2:4));  % Real part
                Q.U(2:4) = q1(1) * q2(2:4) + q2(1) * q1(2:4) + cross(q1(2:4), q2(2:4));  % Imaginary parts
            elseif ~isa(Q1, 'quaternion')
                % Scalar * Quaternion multiplication
                Q = Q2;
                Q.U = Q1 * Q2.U;  % Scalar multiplication
            elseif ~isa(Q2, 'quaternion')
                % Quaternion * Scalar multiplication
                Q = Q1;
                Q.U = Q1.U * Q2;  % Scalar multiplication
            end
        end
        
        % Quaternion addition (overloaded '+')
        function Q = plus(Q1, Q2)
            % Add two quaternions component-wise
            Q = Q1;
            Q.U = Q.U + Q2.U;
        end
        
        % Unary minus (negate the quaternion)
        function Q = uminus(Q)
            % Negate the entire quaternion
            Q.U = -Q.U;
        end
        
        % Quaternion subtraction (overloaded '-')
        function Q = minus(Q1, Q2)
            % Subtract Q2 from Q1 (Q1 - Q2)
            Q = Q1 + (-Q2);  % Reuse addition and negation
        end
        
        % Display the quaternion in symbolic form
        function disp(Q)
            % Display quaternion as a symbolic expression
            expressionForm = Q.expression;
            disp(expressionForm);  % Display the expression
        end
    end
end
