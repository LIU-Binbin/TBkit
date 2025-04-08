function enforce_2D_POSCAR(origin_POSCAR, final_POSCAR, vacuum_length, minimum_bond_length)
    % ENFORCE_2D_POSCAR - Modifies a POSCAR structure to enforce a 2D-like structure
    % by adjusting the lattice vectors and possibly exchanging axes based on a vacuum length.
    %
    % The function reads the original POSCAR, modifies the lattice vectors so that one axis
    % (usually the 'c' axis) is elongated to the specified vacuum length, and outputs the 
    % modified POSCAR structure.
    %
    % Input Arguments:
    %   origin_POSCAR       - The original POSCAR file to be read (string, file path).
    %   final_POSCAR        - The output POSCAR file with enforced 2D structure (string, file path).
    %   vacuum_length       - The desired vacuum length to enforce on the 'c' direction (scalar).
    %   minimum_bond_length - Minimum bond length between atoms to scale the lattice vectors (optional, scalar).
    %
    % Output:
    %   A new POSCAR file is generated with the modified lattice vectors.
    %
    % Example:
    %   enforce_2D_POSCAR('POSCAR_origin', 'POSCAR_final', 10, 1.5);

    % Read the original POSCAR data
    [Rm, sites, Atom_name, Atom_num, ~] = POSCAR_readin(origin_POSCAR);  % Rm: lattice matrix, sites: atom positions
    
    % Get the nearest neighbor distances (used for scaling)
    [~, ~, Rnn, ~] = nn_smart(Rm, sites, [1 1 1], 3);
    minimum_bond = Rnn(1);  % Minimum bond length between atoms
    
    % If the minimum bond length is not provided, use the calculated value
    if nargin < 4
        minimum_bond_length = minimum_bond;
    end
    
    % Scale the lattice matrix (Rm) based on the desired vacuum length
    Rm = Rm * minimum_bond_length / minimum_bond;
    
    % Calculate the variance of each component of the lattice vectors
    Rf_list = [[sites.rc1]', [sites.rc2]', [sites.rc3]'];
    var_a = var(Rf_list(:, 1));  % Variance of the 'a' component
    var_b = var(Rf_list(:, 2));  % Variance of the 'b' component
    var_c = var(Rf_list(:, 3));  % Variance of the 'c' component
    
    % Calculate the lengths of the scaled lattice vectors
    var_a_length = norm(var_a * Rm(1, :));
    var_b_length = norm(var_b * Rm(2, :));
    var_c_length = norm(var_c * Rm(3, :));
    
    Rm_length1 = norm(Rm(1, :));  % Length of the first lattice vector
    Rm_length2 = norm(Rm(2, :));  % Length of the second lattice vector
    Rm_length3 = norm(Rm(3, :));  % Length of the third lattice vector
    
    % Find the smallest variance direction (to choose which axis to scale)
    var_min = min([var_a, var_b, var_c]);
    
    test_label = 0;  % Flag to track if any action is taken
    
    % Determine which direction to scale based on the smallest variance
    switch var_min
        case var_a
            c_dir = 1;  % 'a' direction
            if Rm_length1 >= vacuum_length
                fprintf('Exchange a with c\n');
                test_label = 1;
                C_dir = 1;
            end
        case var_b
            c_dir = 2;  % 'b' direction
            if Rm_length2 >= vacuum_length
                fprintf('Exchange b with c\n');
                test_label = 1;
                C_dir = 2;
            end
        case var_c
            c_dir = 3;  % 'c' direction
            if Rm_length3 >= vacuum_length
                fprintf('No need to change\n');
                test_label = 1;
                C_dir = 3;
            end
        otherwise
            test_label = 2;  % No suitable action
    end
    
    % If no decision was made, use the largest lattice vector
    if test_label == 0
        [~, c_dir] = max([Rm_length1, Rm_length2, Rm_length3]);
        fprintf('Use max length in %d direction\n', c_dir);
    end
    
    % If an action was decided (or default to the largest axis)
    if test_label == 1 || test_label == 0
        switch c_dir
            case 1
                % Scale the 'a' direction to the vacuum length
                if Rm_length1 == max([Rm_length1, Rm_length2, Rm_length3])
                    Rm(1, :) = Rm(1, :) * vacuum_length / Rm_length1;
                    C_dir = 1;
                else
                    fprintf('Program cannot do anything, take your own action!\n');
                    return;
                end
            case 2
                % Scale the 'b' direction to the vacuum length
                if Rm_length2 == max([Rm_length1, Rm_length2, Rm_length3])
                    Rm(2, :) = Rm(2, :) * vacuum_length / Rm_length2;
                    C_dir = 2;
                else
                    fprintf('Program cannot do anything, take your own action!\n');
                    return;
                end
            case 3
                % Scale the 'c' direction to the vacuum length
                if Rm_length3 == max([Rm_length1, Rm_length2, Rm_length3])
                    Rm(3, :) = Rm(3, :) * vacuum_length / Rm_length3;
                    C_dir = 3;
                else
                    fprintf('Program cannot do anything, take your own action!\n');
                    return;
                end
        end
    end
    
    % Apply the correct permutation matrix based on the chosen axis for scaling
    switch C_dir
        case 3
            % No need to exchange axes, just generate the final POSCAR
            POSCAR_gen(Rm, sites, Atom_name, Atom_num, final_POSCAR);
        case 2
            % Exchange the 'b' and 'c' directions
            P_matrix = [1 0 0; 0 0 1; 0 1 0];
            Rm = P_matrix * Rm * P_matrix.';  % Apply permutation to the lattice matrix
            Rf_list = Rf_list * P_matrix;    % Apply the permutation to the atomic positions
            for i = 1:length(sites)
                sites(i).rc1 = Rf_list(i, 1);
                sites(i).rc2 = Rf_list(i, 2);
                sites(i).rc3 = Rf_list(i, 3);
            end
            POSCAR_gen(Rm, sites, Atom_name, Atom_num, final_POSCAR);
        case 1
            % Exchange the 'a' and 'c' directions
            P_matrix = [0 0 1; 0 1 0; 1 0 0];
            Rm = P_matrix * Rm * P_matrix.';  % Apply permutation to the lattice matrix
            Rf_list = Rf_list * P_matrix;    % Apply the permutation to the atomic positions
            for i = 1:length(sites)
                sites(i).rc1 = Rf_list(i, 1);
                sites(i).rc2 = Rf_list(i, 2);
                sites(i).rc3 = Rf_list(i, 3);
            end
            POSCAR_gen(Rm, sites, Atom_name, Atom_num, final_POSCAR);  % Generate the final POSCAR file
    end
end
