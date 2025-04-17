function [Rm, sites] = supercell(Ns, Rm, sites, Atom_name, Atom_num, findir, filename)
    %% Supercell Generation Function
    % This function generates a supercell based on the input parameters. It computes the 
    % new atomic positions for the supercell and writes the results to a POSCAR file.
    % 
    % Usage:
    %   [Rm, sites] = supercell(Ns, Rm, sites, Atom_name, Atom_num, findir, filename)
    %
    % Inputs:
    %   Ns - Supercell matrix (3x3), defines the number of unit cells in each direction
    %   Rm - Lattice matrix (3x3) for the unit cell
    %   sites - Struct array containing the atomic positions in reduced coordinates
    %   Atom_name - Cell array of atomic element names
    %   Atom_num - Array containing the number of atoms for each element
    %   findir - Optional, directional shift vector for the supercell (default: [0, 0, 0])
    %   filename - Optional, name of the output POSCAR file (default: 'POSCAR_super')
    %
    % Outputs:
    %   Rm - Updated lattice matrix for the supercell
    %   sites - Struct array of updated atomic positions for the supercell
    
    %% Default Arguments Handling
    if nargin < 2
        [Rm, sites, Atom_name, Atom_num, elements, a_crystal_constance] = POSCAR_read;
        findir = [0, 0, 0];  % Default shift vector
        filename = 'POSCAR_super';  % Default output file name
    end
    if nargin < 6
        findir = [0, 0, 0];  % Default shift vector
        filename = 'POSCAR_super';  % Default output file name
    end
    if nargin < 7
        filename = 'POSCAR_super';  % Default output file name
    end
    
    %% Calculate the Volume of the Supercell
    V_Ns = dot(Ns(1,:), cross(Ns(2,:), Ns(3,:)));  % Volume of the supercell
    if V_Ns < 0
        Ns(3,:) = -Ns(3,:);  % Ensure positive volume
    end
    V_Ns = abs(V_Ns);
    
    if V_Ns < 1
        error('The supercell volume is too small. Check Ns input.');
    end
    
    %% Adjust Atomic Positions (Pre-processing)
    % Update atomic positions using the plusrc function (wrap coordinates within [0, 1))
    [~, nsites] = size(sites);
    for i = 1:nsites
        sites(i).rc1 = plusrc(sites(i).rc1);
        sites(i).rc2 = plusrc(sites(i).rc2);
        sites(i).rc3 = plusrc(sites(i).rc3);
    end
    
    %% Generate Supercell Vectors and Atomic Positions
    max_R = max(abs(Ns)) * 3;  % Maximum range for generating supercell vectors
    sc_cands = generateSupercellCandidates(max_R);  % Generate candidate translation vectors
    
    % Initialize supercell vectors
    sc_vec = [];
    eps_shift = sqrt(2.0) * 1.0E-8;  % Small epsilon to avoid double counting
    
    % Check all candidate vectors to see if they are inside the supercell
    for ivec = 1:length(sc_cands)
        tmp_red = to_red_sc(sc_cands(ivec,:), Ns);  % Convert to reduced coordinates
        inside = all(tmp_red >= -eps_shift & tmp_red < 1 - eps_shift);  % Check if inside the unit cell
        
        if inside
            sc_vec = [sc_vec; sc_cands(ivec,:)];  % Add valid vector to supercell
        end
    end
    
    % Ensure correct number of supercell vectors
    [num_sc, ~] = size(sc_vec);
    if round(abs(det(Ns))) ~= num_sc
        error('Supercell generation failed! Wrong number of supercell vectors found.');
    end
    
    %% Compute New Atomic Positions for Supercell
    countnum = 0;
    count = 0;
    for elementseq = 1:length(Atom_num)
        for n = 1:Atom_num(elementseq)
            countnum = countnum + 1;
            Rpc = [sites(countnum).rc1, sites(countnum).rc2, sites(countnum).rc3];  % Original position
            % Loop over each supercell vector
            for icur_sc_vec = 1:num_sc
                count = count + 1;
                cur_sc_vec = sc_vec(icur_sc_vec, :);  % Current supercell translation vector
                Rsc = to_red_sc(Rpc + cur_sc_vec, Ns);  % Compute reduced coordinates for new position
                Rsc = round(Rsc * 10^8) / 10^8;  % Round to avoid floating point errors
                % Ensure atomic coordinates stay within [0, 1)
                Rsc = mod(Rsc, 1);
                sites_s(count).rc1 = Rsc(1);
                sites_s(count).rc2 = Rsc(2);
                sites_s(count).rc3 = Rsc(3);
                sites_s(count).name = Atom_name(elementseq);
            end
        end
    end
    sites = sites_s;  % Update atomic positions
    Rm = Ns * Rm;
    %% Apply Directional Shift (if findir is provided)
    if any(findir)
        Rm_s_fin = applyFindirShift(Rm, findir);  % Apply shift to lattice vectors
        sites = shiftAtomicPositions(sites, Rm_s_fin);  % Apply shift to atomic positions
        Rm = Rm_s_fin;  % Update lattice matrix
    end
    
    %% Generate POSCAR File
    Atom_num_s = round(Atom_num * V_Ns);  % Update atom counts for supercell
    POSCAR_gen(Rm, sites, Atom_name, Atom_num_s, filename);  % Generate POSCAR file
    
end

%% Helper Functions

% Function to wrap coordinates within [0, 1)
function rc_plus = plusrc(rc)
    if rc < 0
        rc_plus = rc + 1;  % Wrap negative values to [0, 1)
    elseif rc >= 1
        rc_plus = mod(rc, 1);  % Wrap values >= 1 to [0, 1)
    else
        rc_plus = rc;  % No change needed
    end
end

% Function to generate supercell translation candidate vectors
function sc_cands = generateSupercellCandidates(max_R)
    sc_cands = [];
    for i = -max_R(1):max_R(1)
        for j = -max_R(2):max_R(2)
            for k = -max_R(3):max_R(3)
                sc_cands = [sc_cands; [i, j, k]];  % Store translation vectors
            end
        end
    end
end

% Function to apply directional shift (findir) to supercell
function Rm_s_fin = applyFindirShift(Rm, findir)
    % Adjust the lattice vectors based on findir
    Rmlength1 = abs(norm(Rm(1,:)));
    Rmlength2 = abs(norm(Rm(2,:)));
    Rmlength3 = abs(norm(Rm(3,:)));
    
    Rm_s_fin_add = [10 * Rm(1,:) * findir(1) / Rmlength1;
                    10 * Rm(2,:) * findir(2) / Rmlength2;
                    10 * Rm(3,:) * findir(3) / Rmlength3];
    Rm_s_fin = Rm + Rm_s_fin_add;
end

% Function to shift atomic positions based on new lattice vectors
function sites = shiftAtomicPositions(sites, Rm_s_fin)
    % Apply the shift to the atomic positions
    for i = 1:length(sites)
        Rr = [sites(i).rc1, sites(i).rc2, sites(i).rc3] * Rm_s_fin;  % Cartesian coordinates
        Rc_s_fin = Rr / Rm_s_fin;  % Convert back to reduced coordinates
        sites(i).rc1 = Rc_s_fin(1);
        sites(i).rc2 = Rc_s_fin(2);
        sites(i).rc3 = Rc_s_fin(3);
    end
end
