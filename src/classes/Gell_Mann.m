classdef Gell_Mann < pauli_matrix
    % Gell_Mann class for handling Gell-Mann matrices.
    % This class extends the pauli_matrix class and constructs 
    % Gell-Mann matrices, commonly used in quantum mechanics, 
    % especially in the context of SU(3) symmetries.

    properties
        % No additional properties are required, as the class inherits from pauli_matrix
    end
    
    methods
        function GM = Gell_Mann(label1, label2, options)
            % Constructor for the Gell_Mann class.
            % This creates Gell-Mann matrices based on the given labels and options.
            
            arguments
                label1 double = 0;  % First label for selecting the matrix
                label2 double = 0;  % Second label (used in commutation)
                options.rep char = 'norm';  % Representation style, default is 'norm'
            end
            
            % Call the parent constructor (pauli_matrix)
            GM = GM@pauli_matrix(0);  
            paulis = pauli_matrix();  % Create an instance of the Pauli matrices
            pauli_mat = zeros(3);     % Initialize a 3x3 matrix for the Gell-Mann matrix
            
            % Handle the case where only label1 is provided
            if nargin == 1
                switch options.rep
                    case 'norm'  % Normal representation (SU(3) specific)
                        % Define Gell-Mann matrices based on the label1
                        pauli_mat = defineGellMannMatrix(label1, paulis);
                    otherwise
                        % Alternative representation (unspecified)
                        pauli_mat = defineGellMannMatrix(label1, paulis);
                end
                GM.mat = pauli_mat;  % Assign the generated matrix to the class property
            elseif nargin == 0
                % Generate a full set of 8 Gell-Mann matrices if no arguments are provided
                count = 0;
                for i = 1:8
                    count = count + 1;
                    GM(count) = Gell_Mann(i, 'rep', options.rep);  % Generate and store each matrix
                end
            else
                % If two labels are provided, calculate the commutation of two matrices
                GM = gamma_matrix(label1, 'rep', options.rep);
                GM = GM.commute(Gell_Mann(label2, 'rep', options.rep));
                GM.mat = GM.mat / (2 * 1i);  % Normalize the result by dividing by 2i
            end
        end
    end
    
    methods (Static)
        % Add any static methods if needed (currently empty)
    end
end

function pauli_mat = defineGellMannMatrix(label, paulis)
    % Helper function to define a Gell-Mann matrix based on the given label.
    % This function simplifies the matrix selection code and improves readability.
    
    pauli_mat = zeros(3);  % Initialize a 3x3 zero matrix

    switch label
        case 0
            pauli_mat = eye(3);  % Identity matrix for label 0
        case 1
            pauli_mat(1:2, 1:2) = paulis(2).mat;  % Pauli matrix 1 in top-left corner
        case 2
            pauli_mat(1:2, 1:2) = paulis(3).mat;  % Pauli matrix 2 in top-left corner
        case 3
            pauli_mat(1:2, 1:2) = paulis(4).mat;  % Pauli matrix 3 in top-left corner
        case 4
            pauli_mat([1, 3], [1, 3]) = paulis(2).mat;  % Pauli matrix 1 in 1st and 3rd rows/columns
        case 5
            pauli_mat([1, 3], [1, 3]) = paulis(3).mat;  % Pauli matrix 2 in 1st and 3rd rows/columns
        case 6
            pauli_mat(2:3, 2:3) = paulis(2).mat;  % Pauli matrix 1 in bottom-right corner
        case 7
            pauli_mat(2:3, 2:3) = paulis(3).mat;  % Pauli matrix 2 in bottom-right corner
        case 8
            pauli_mat = diag([1, 1, -2] / sqrt(3));  % Special case for Gell-Mann matrix 8
    end
end
