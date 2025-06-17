classdef gamma_matrix < pauli_matrix
    % Gamma Matrix class for handling gamma matrices.
    % This class is used to generate gamma matrices in various representations 
    % (such as 'ZYX', 'BJY', 'Dirac', etc.) based on given labels.
    
    properties
        % No additional properties are required, as the class inherits from pauli_matrix
    end
    
    methods
        function GM = gamma_matrix(label1, label2, options)
            % Constructor for the gamma_matrix class.
            % This creates gamma matrices based on the given labels and options.
            
            arguments
                label1 double = 0;   % First label for selecting the matrix
                label2 double = 0;   % Second label (used in commutation)
                options.rep char {mustBeMember(options.rep,{'ZYX','BJY','ZSC','ZSC_2004','KM','Dirac','Weyl'})}= 'ZYX';  % Representation style, default is 'ZYX'
            end
            
            % Call the parent constructor (pauli_matrix)
            GM = GM@pauli_matrix(0);
            paulis = pauli_matrix();  % Create an instance of the Pauli matrices
            
            % Case for single label input (nargin == 1)
            if nargin == 1
                % Switch to select different representations
                switch options.rep
                    case 'ZYX'
                        pauli = selectMatrix(label1, paulis);
                    case 'BJY'
                        pauli = selectMatrix(label1, paulis);
                    case 'ZSC'
                        pauli = selectMatrix(label1, paulis);
                    case 'ZSC_2004'
                        pauli = selectMatrix(label1, paulis);
                    case 'KM'
                        pauli = selectMatrix(label1, paulis);
                    case 'Dirac'
                        pauli = selectMatrix(label1, paulis);
                    case 'Weyl'
                        pauli = selectMatrix(label1, paulis);
                    otherwise
                        pauli = selectMatrix(label1, paulis);
                end
                GM.mat = pauli.mat;  % Assign the selected matrix to the class property
            elseif nargin == 0
                % Generate a full set of gamma matrices if no arguments are provided
                GM = generateFullSet(options.rep);
            else
                % Calculate the commutation of two gamma matrices
                GM = gamma_matrix(label1, 'rep', options.rep);
                GM = GM.commute(gamma_matrix(label2, 'rep', options.rep));
                GM.mat = GM.mat / (2 * 1i);  % Normalize the result by dividing by 2i
            end
        end
    end
    
    methods (Static)
        % Function to generate the inverse of the matrix S
        function Smat_inv = S()
            Smat = zeros(16);
            GM_L = sym(gamma_matrix);
            for i = 1:16
                tmp_mat = GM_L(:,:,i);
                tmp_mat_r = real(tmp_mat);
                tmp_mat_i = imag(tmp_mat);
                Smat(i, 1:4) = diag(tmp_mat_r);
                Smat(i, 5:7) = diag(tmp_mat_r, 1);
                Smat(i, 8:9) = diag(tmp_mat_r, 2);
                Smat(i, 10)  = diag(tmp_mat_r, 3);
                Smat(i, 11:13) = diag(tmp_mat_i, 1);
                Smat(i, 14:15) = diag(tmp_mat_i, 2);
                Smat(i, 16)  = diag(tmp_mat_i, 3);
            end
            Smat_inv = inv(Smat);  % Return the inverse of Smat
        end
        
        % Function to generate a symbolic list of gamma matrices L
        function Gamma_L = L()
            count = 0;
            Gamma_L = sym(zeros(1, 16));
            for i = 0:5
                count = count + 1;
                Gamma_L(count) = sym(['Gamma_', num2str(i)]);
            end
            
            for i = 2:5
                for j = 1:i-1
                    count = count + 1;
                    Gamma_L(count) = sym(['Gamma_', num2str(j), '_', num2str(i)]);
                end
            end
        end
        
        % Function to generate a symbolic set of Pauli matrices
        function Pauli_L = pauli_L()
            Pauli_L = sym(zeros(1, 16));
            syms sigma_0 sigma_x sigma_z tau_0 tau_x tau_z real;
            syms sigma_y tau_y;
            
            % Define the Pauli matrices using symbolic variables
            Pauli_L(1) = sigma_0 * tau_0;
            Pauli_L(2) = sigma_0 * tau_x;
            Pauli_L(3) = sigma_y * tau_y;
            Pauli_L(4) = sigma_0 * tau_z;
            Pauli_L(5) = sigma_x * tau_y;
            Pauli_L(6) = sigma_z * tau_y;
            Pauli_L(7) = -sigma_y * tau_z;
            Pauli_L(8) = sigma_0 * tau_y;
            Pauli_L(9) = -sigma_y * tau_x;
            Pauli_L(10) = -sigma_x * tau_z;
            Pauli_L(11) = sigma_z * tau_0;
            Pauli_L(12) = sigma_x * tau_x;
            Pauli_L(13) = -sigma_z * tau_z;
            Pauli_L(14) = -sigma_x * tau_0;
            Pauli_L(15) = sigma_z * tau_x;
            Pauli_L(16) = sigma_y * tau_0;
        end
    end
end

% Helper function for matrix selection based on label and Pauli matrices
function pauli = selectMatrix(label, paulis)
    switch label
        case 0
            pauli = paulis(1) * paulis(1);
        case 1
            pauli = paulis(1) * paulis(2);
        case 2
            pauli = paulis(3) * paulis(3);
        case 3
            pauli = paulis(1) * paulis(4);
        case 4
            pauli = paulis(2) * paulis(3);
        case 5
            pauli = paulis(4) * paulis(3);
    end
end

% Helper function to generate the full set of gamma matrices based on representation
function GM = generateFullSet(rep)
    count = 0;
    GM = gamma_matrix(0, 'rep', rep);
    
    switch rep
        case {'Dirac', 'Weyl'}
            for i = [0, 1, 2, 3, 5]
                for j = [0, 1, 2, 3, 5]
                    if i < j
                        count = count + 1;
                        GM(count) = gamma_matrix(i, j, 'rep', rep);
                    end
                end
            end
        otherwise
            for i = 2:5
                for j = 1:i-1
                    count = count + 1;
                    GM(count) = gamma_matrix(i, j, 'rep', rep);
                end
            end
    end
end
