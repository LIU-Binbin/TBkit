function H_hr = cut_piece(H_hr, repeatnum, fin_dir, glue_edges, vacuum_mode)
% CUT_PIECE Constructs a supercell along a specified direction with optional vacuum layer and edge conditions
%   This function modifies an HR object to create a supercell by repeating the unit cell along a given direction.
%   It handles different Hamiltonian matrix formats ('sparse', 'mat', 'list') and supports vacuum layer addition
%   as well as periodic (glued) or non-periodic edge conditions.
%
% INPUTS:
%   H_hr       - HR object containing Hamiltonian and structural information
%   repeatnum  - Number of repetitions along the specified direction (default: 10, must be positive integer)
%   fin_dir    - Direction of repetition (1: x, 2: y, 3: z, 4: custom; default: 3)
%   glue_edges - Logical flag for periodic boundary conditions (glue edges, default: false)
%   vacuum_mode - Logical flag to add vacuum layer along the repetition direction (default: false)
%
% OUTPUTS:
%   H_hr       - Modified HR object with updated supercell Hamiltonian and structural properties
%
% NOTES:
%   - repeatnum must be ≥ 1; glue_edges is invalid when repeatnum = 1
%   - vacuum_mode adds a 10x unit cell length vacuum layer along fin_dir
%   - Handles different Hamiltonian storage types (sparse matrices, full matrices, list format)
%   - Validates input arguments and throws errors for invalid configurations
%   - Updates orbital positions, quantum numbers, and element labels for the supercell

arguments
    H_hr HR;
    repeatnum double{mustBeInteger, mustBePositive} = 10;
    fin_dir double{mustBeMember(fin_dir, [1, 2, 3, 4])} = 3;
    glue_edges logical = false;
    vacuum_mode logical = false;
end

% Validate input combinations
if repeatnum == 1 && glue_edges
    error("\n\nCan't have num==1 and glueing of the edges!");
end

% Construct supercell expansion matrix
Ns = eye(3);
Ns(fin_dir, :) = Ns(fin_dir, :) * repeatnum;

% Generate supercell orbital information
[sc_orb, ~, sc_elementL, sc_quantumL] = H_hr.supercell_orb(Ns);
sc_orb = double(sc_orb);

% Handle vacuum layer addition
if vacuum_mode
    fin_dir_list = zeros(1, 3);
    fin_dir_list(fin_dir) = 1; % Mark vacuum direction
    
    % Calculate vacuum layer lattice vector addition
    Rmlength = norm(H_hr.Rm(fin_dir, :));
    Rm_s_fin_add = 10 * H_hr.Rm(fin_dir, :) / Rmlength .* fin_dir_list;
    Rm_s_fin = H_hr.Rm + Rm_s_fin_add;
    
    % Adjust orbital positions to new supercell coordinates
    Rr_s_fin_add = [0.5, 0.5, 0.5] * Rm_s_fin_add; % Center vacuum layer
    for i = 1:size(sc_orb, 1)
        Rr_orb = sc_orb(i, :) * H_hr.Rm;
        Rr_s_fin = Rr_orb + Rr_s_fin_add;
        sc_orb(i, :) = Rr_s_fin / Rm_s_fin; % Convert to fractional coordinates
    end
end

% Generate supercell structure from POSCAR if sites are not present
if isempty(H_hr.sites)
    if exist('POSCAR', 'file')
        [Rm_tmp, sites_tmp, Atom_name_tmp, Atom_num_tmp, ~] = HR.POSCAR_read('POSCAR', 'vasp');
        try
            H_hr = H_hr.supercell(Ns, 'POSCAR_super_fin', Rm_tmp, sites_tmp, Atom_name_tmp, Atom_num_tmp, fin_dir_list);
        catch ME
            warning('Supercell generation from POSCAR failed: %s', ME.message);
        end
end

% Hamiltonian transformation based on storage type
WANNUM = H_hr.WAN_NUM;
OUT_WAN_NUM = WANNUM * repeatnum;

switch H_hr.Type
    case 'sparse'
        % Initialize sparse Hamiltonian blocks
        OUT_HnumL = cell(H_hr.NRPTS, 1);
        for i = 1:H_hr.NRPTS
            OUT_HnumL{i} = sparse(OUT_WAN_NUM, OUT_WAN_NUM);
        end
        
        pb = TBkit_tool_outer.CmdLineProgressBar('NRPT : ');
        for iN = 1:H_hr.NRPTS
            ind_R = H_hr.vectorL(iN, :);
            jump_fin = ind_R(fin_dir);
            
            % Process each hopping term
            [ilist, jlist, amplist] = find(H_hr.HnumL{iN});
            for ih = 1:length(amplist)
                for icur_sc_vec = 1:repeatnum
                    hi = ilist(ih) + (icur_sc_vec - 1) * WANNUM;
                    hj = jlist(ih) + (icur_sc_vec + jump_fin - 1) * WANNUM;
                    
                    % Apply boundary conditions
                    if glue_edges
                        hj = mod(hj, OUT_WAN_NUM);
                        hj = hj == 0 ? OUT_WAN_NUM : hj;
                    else
                        if hj < 1 || hj > OUT_WAN_NUM, continue; end
                    end
                    
                    OUT_HnumL{iN}(hi, hj) = amplist(ih);
                end
            end
            pb.print(iN, H_hr.NRPTS);
        end
        pb.delete();
        
        H_hr.HnumL = OUT_HnumL;
        
    case 'mat'
        % Handle full matrix representation
        OUT_HnumL = zeros(OUT_WAN_NUM, OUT_WAN_NUM, H_hr.NRPTS);
        vectorList = double(H_hr.vectorL);
        
        pb = TBkit_tool_outer.CmdLineProgressBar('NRPT : ');
        for ih = 1:H_hr.NRPTS
            ind_R = vectorList(ih, :);
            jump_fin = ind_R(fin_dir);
            tmpHnum = H_hr.HnumL(:, :, ih);
            
            % Map original indices to supercell indices
            [hi, hj] = find(tmpHnum);
            for icur_sc_vec = 1:repeatnum
                hi_sc = hi + (icur_sc_vec - 1) * WANNUM;
                hj_sc = hj + (icur_sc_vec + jump_fin - 1) * WANNUM;
                
                % Apply boundary conditions and accumulate contributions
                if glue_edges
                    hj_sc = mod(hj_sc, OUT_WAN_NUM);
                    hj_sc(hj_sc == 0) = OUT_WAN_NUM;
                else
                    valid = (hj_sc > 0) & (hj_sc <= OUT_WAN_NUM);
                    hi_sc = hi_sc(valid);
                    hj_sc = hj_sc(valid);
                end
                
                OUT_HnumL(hi_sc, hj_sc, ih) = tmpHnum(sub2ind([WANNUM, WANNUM], hi(valid), hj(valid)));
            end
            pb.print(ih, H_hr.NRPTS);
        end
        pb.delete();
        
        H_hr.HnumL = OUT_HnumL;
        
    case 'list'
        % Process list-based Hamiltonian representation
        OUT_vectorList = [];
        OUT_HnumL = [];
        
        pb = TBkit_tool_outer.CmdLineProgressBar('H:NRPT : ');
        for ih = 1:H_hr.NRPTS
            ind_R = H_hr.vectorL(ih, 1:H_hr.Dim);
            jump_fin = ind_R(fin_dir);
            i = H_hr.vectorL(ih, H_hr.Dim + 1);
            j = H_hr.vectorL(ih, H_hr.Dim + 2);
            amp = H_hr.HnumL(ih);
            
            for icur_sc_vec = 1:repeatnum
                hi = i + (icur_sc_vec - 1) * WANNUM;
                hj = j + (icur_sc_vec + jump_fin - 1) * WANNUM;
                
                if glue_edges
                    hj = mod(hj, OUT_WAN_NUM);
                    hj = hj == 0 ? OUT_WAN_NUM : hj;
                else
                    if hj < 1 || hj > OUT_WAN_NUM, continue; end
                end
                
                OUT_vectorList = [OUT_vectorList; [ind_R, hi, hj]];
                OUT_HnumL = [OUT_HnumL; amp];
            end
            pb.print(ih, H_hr.NRPTS);
        end
        pb.delete();
        
        H_hr.vectorL = int8(OUT_vectorList);
        H_hr.HnumL = OUT_HnumL;
end

% Finalize supercell properties
H_hr.orbL = sc_orb;
H_hr.quantumL = sc_quantumL;
H_hr.elementL = sc_elementL;
H_hr.num = true;
H_hr.coe = false;
end