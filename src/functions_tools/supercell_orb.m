function orbital_out = supercell_orb(orbital_init, fin_dir, vacuum_mode)
    %% Function to generate supercell orbital positions
    % This function computes the final orbital positions after constructing 
    % the supercell, depending on the vacuum_mode.
    % 
    % Usage:
    %   orbital_out = supercell_orb(orbital_init, fin_dir, vacuum_mode)
    %
    % Inputs:
    %   orbital_init - Initial orbital positions (Nx3 matrix)
    %   fin_dir - Final direction to expand the supercell (1x3 vector)
    %   vacuum_mode - Mode to handle vacuum space in the supercell (0 or 1)
    %
    % Output:
    %   orbital_out - Final orbital positions after supercell construction
    
    switch vacuum_mode
        case 0
            % Case 0: Supercell construction without vacuum
            orbital_out = orbital_init;
            Ns = [1 0 0; 0 1 0; 0 0 1]; 
            Ns = Ns .* fin_dir;  % Scale the supercell matrix by fin_dir
            
            fin_dir_list = [0 0 0];  % No vacuum direction
            
            % Read POSCAR and generate supercell
            [Rm, sites, Atom_name, Atom_num] = POSCAR_read();
            [~, ~] = supercell(Ns, Rm, sites, Atom_name, Atom_num, fin_dir_list, 'POSCAR_super_fin');
            
            % Loop through each dimension for finite directions and modify orbital positions
            for i = 1:3
                Nslab = fin_dir(i);
                if Nslab == 0
                    Nslab = 1;  % Ensure non-zero slab count
                end
                
                % Initialize counter and orbital storage
                count = 0;
                WAN_NUM = size(orbital_out, 1);
                fin_orb = zeros(WAN_NUM * Nslab, 3);  % Initialize final orbital positions
                
                % Loop over all finite slabs and orbital positions
                for inum = 1:Nslab
                    for j = 1:WAN_NUM
                        count = count + 1;
                        
                        % Copy the j-th orbital and adjust its coordinate in the i-th direction
                        orb_tmp = orbital_out(j, :);
                        orb_tmp(i) = (orb_tmp(i) + double(inum - 1)) / Nslab;  % Adjust fractional coordinate
                        fin_orb(count, :) = orb_tmp;  % Store the modified orbital
                    end
                end
                orbital_out = fin_orb;  % Update orbital positions
            end
            
        case 1
            % Case 1: Supercell construction with vacuum handling
            orbital_out = orbital_init;
            
            % Loop through each direction and expand the supercell
            for i = 1:3
                Nslab = fin_dir(i);
                if Nslab == 0
                    Nslab = 1;  % Ensure non-zero slab count
                end
                
                count = 0;
                WAN_NUM = size(orbital_out, 1);
                fin_orb = zeros(WAN_NUM * Nslab, 3);  % Initialize final orbital positions
                
                % Loop over all finite slabs and orbital positions
                for inum = 1:Nslab
                    for j = 1:WAN_NUM
                        count = count + 1;
                        
                        % Copy the j-th orbital and adjust its coordinate in the i-th direction
                        orb_tmp = orbital_out(j, :);
                        orb_tmp(i) = (orb_tmp(i) + double(inum - 1)) / Nslab;  % Adjust fractional coordinate
                        fin_orb(count, :) = orb_tmp;  % Store the modified orbital
                    end
                end
                orbital_out = fin_orb;  % Update orbital positions
            end
            
            % Generate the supercell (POSCAR generation)
            Ns = [1 0 0; 0 1 0; 0 0 1]; 
            Ns = Ns .* fin_dir;
            fin_dir_list = double(fin_dir > 1);  % Determine the direction where vacuum will be
            
            % Read POSCAR and generate supercell
            [Rm, sites, Atom_name, Atom_num] = POSCAR_read();
            [~, ~] = supercell(Ns, Rm, sites, Atom_name, Atom_num, fin_dir_list, 'POSCAR_super_fin');
            
            % Modify lattice matrix based on the supercell
            Rm = Ns * Rm;
            Rmlength1 = norm(Rm(1, :));
            Rmlength2 = norm(Rm(2, :));
            Rmlength3 = norm(Rm(3, :));
            
            % Calculate lattice adjustments based on fin_dir (vacuum handling)
            Rm_s_fin_add = [10 * Rm(1, :) * fin_dir_list(1) / Rmlength1;
                            10 * Rm(2, :) * fin_dir_list(2) / Rmlength2;
                            10 * Rm(3, :) * fin_dir_list(3) / Rmlength3];
            Rm_s_fin = Rm + Rm_s_fin_add;  % Adjusted lattice matrix
            Rc_s_fin_add = [1 / 2, 1 / 2, 1 / 2];  % Shift for vacuum direction
            Rr_s_fin_add = Rc_s_fin_add * Rm_s_fin_add;  % Adjusted fractional position for vacuum
            
            % Loop through all orbital positions and apply the supercell transformation
            [nfinorb, ~] = size(fin_orb);
            for i = 1:nfinorb
                Rr_orb = fin_orb(i, :) * Rm;  % Cartesian coordinates of orbital
                Rr_s_fin = Rr_orb + Rr_s_fin_add;  % Adjusted coordinates
                Rc_s_fin = Rr_s_fin / Rm_s_fin;  % Convert back to reduced coordinates
                fin_orb(i, :) = Rc_s_fin;  % Store the final orbital position
            end
            
            orbital_out = fin_orb;  % Update orbital positions
    end
end
