function H_hr = cut_piece(H_hr, repeatnum, fin_dir, glue_edges, vacuum_mode, Accuracy)
%CUT_PIECE Construct lower-dimensional model from periodic Hamiltonian via dimensional reduction
%
%   H_hr = CUT_PIECE(H_hr, repeatnum, fin_dir, glue_edges, vacuum_mode, Accuracy)
%   Generates a (d-1)-dimensional tight-binding model by repeating a d-dimensional
%   Hamiltonian along specified lattice direction. Supports edge gluing and vacuum layer insertion.
%
%   Features:
%   - Handles sparse/matrix/list-type Hamiltonian formats
%   - Maintains orbital positions and quantum numbers
%   - Supports periodic boundary conditions (edge gluing)
%   - Optional vacuum layer insertion for surface systems
%
%   Input Parameters:
%     H_hr         (HR)      : Original Hamiltonian object with periodic boundary conditions
%     repeatnum    (integer) : Number of unit cell repetitions along cutting direction (default=10)
%     fin_dir      (1-3)     : Lattice vector index for dimensional reduction (1=x,2=y,3=z) (default=3)
%     glue_edges   (logical) : Enable periodic boundary along cutting direction (default=false)
%     vacuum_mode  (logical) : Insert vacuum layer in cutting direction (default=false)
%     Accuracy     (double)  : Numerical tolerance for supercell construction (default=1e-6)
%
%   Output:
%     H_hr         (HR)      : Modified Hamiltonian object with (d-1)-dimensional structure
%
%   Notes:
%     - When vacuum_mode=true, the cutting direction lattice vector is expanded by 10x
%     - Edge gluing (glue_edges=true) implements periodic boundary conditions
%     - For fin_dir=4, user must provide custom lattice vectors via H_hr.Rm
%
%   Example:
%     % Create 2D slab from 3D Hamiltonian with 5 layers
%     H_2D = cut_piece(H_3D, 5, 3, false, true);

% ---- Argument Validation ----
arguments
    H_hr HR
    repeatnum double{mustBeInteger} = 10
    fin_dir double{mustBeMember(fin_dir,[1,2,3,4])} = 3
    glue_edges logical = false
    vacuum_mode logical = false
    Accuracy = 1e-6
end

% ---- Supercell Construction ----
% Generate supercell by repeating along specified direction
Ns = [1 0 0; 0 1 0; 0 0 1]; % Default repetition factors
Ns(fin_dir,:) = Ns(fin_dir,:) * repeatnum; % Apply repetition along cutting direction

% Generate supercell orbitals and lattice vectors
[sc_orb, sc_vec, sc_elementL, sc_quantumL] = H_hr.supercell_orb(Ns, Accuracy);
sc_orb = double(sc_orb); % Ensure numerical precision

% ---- POSCAR Handling ----
% Update crystal structure information if vacuum layer requested
fin_dir_list = [0 0 0]; % Direction flags for vacuum expansion
if vacuum_mode
    fin_dir_list(fin_dir) = 1; % Mark direction for vacuum expansion
end

% Transfer crystal structure data to new supercell
if isempty(H_hr.sites)
    if exist('POSCAR','file') % Load from VASP file if available
        [Rm_tmp, sites_tmp, Atom_name_tmp, Atom_num_tmp, ~] = POSCAR_read('POSCAR','vasp');
        H_hr = H_hr.supercell(Ns, 'POSCAR_super_fin', Rm_tmp, sites_tmp, ...
            Atom_name_tmp, Atom_num_tmp, fin_dir_list);
    else
        warning('Crystal structure information not available');
    end
else % Use existing crystal data
    try
        H_hr = H_hr.supercell(Ns, 'POSCAR_super_fin', H_hr.Rm, H_hr.sites, ...
            H_hr.Atom_name, H_hr.Atom_num, fin_dir_list);
    catch ME
        warning('Supercell construction failed: %s', ME.message);
    end
end

% ---- Orbital Consistency Check ----
if size(H_hr.orbL,1) ~= H_hr.WAN_NUM
    error("Orbital initialization mismatch: orbitalL dimension != WAN_NUM");
end

% ---- Input Validation ----
if repeatnum < 1
    error("Repetition number must be positive");
end
if repeatnum == 1 && glue_edges
    error("Single repetition incompatible with edge gluing");
end

% ---- Vacuum Layer Processing ----
if vacuum_mode
    % Expand lattice vectors in cutting direction
    Rm = H_hr.Rm;
    expansion_factors = 10*fin_dir_list./vecnorm(Rm,2,2)';
    Rm_s_fin_add = diag(expansion_factors) * Rm;
    Rm_s_fin = Rm + Rm_s_fin_add;
    
    % Adjust orbital positions for expanded lattice
    Rc_s_fin_add = [0.5, 0.5, 0.5]; % Center shift for vacuum layer
    Rr_s_fin_add = Rc_s_fin_add * Rm_s_fin_add;
    
    % Convert orbital coordinates to new lattice basis
    for i = 1:size(sc_orb,1)
        Rr_orb = sc_orb(i,:) * Rm;
        Rr_s_fin = Rr_orb + Rr_s_fin_add;
        sc_orb(i,:) = Rr_s_fin / Rm_s_fin;
    end
end

% ---- Hamiltonian Reshaping ----
OUT_WAN_NUM = H_hr.WAN_NUM * repeatnum; % New system size
vectorList = double(H_hr.vectorL); % Hopping vectors
NRPT_seq = false(H_hr.NRPTS,1); % Hopping selection mask

% Process different Hamiltonian storage formats
switch H_hr.Type
    case 'sparse'
        % Sparse matrix format processing
        H_hr = process_sparse(H_hr, OUT_WAN_NUM, repeatnum, fin_dir, glue_edges, vectorList);
        
    case 'mat'
        % Full matrix format processing
        H_hr = process_matrix(H_hr, OUT_WAN_NUM, repeatnum, fin_dir, glue_edges, vectorList, NRPT_seq);
        
    case 'list'
        % List format processing
        H_hr = process_list(H_hr, OUT_WAN_NUM, repeatnum, fin_dir, glue_edges, vectorList);
        
    otherwise
        error('Unsupported Hamiltonian type: %s', H_hr.Type);
end

% ---- Finalize HR Object ----
H_hr.orbL = sc_orb; % Update orbital positions
H_hr.quantumL = sc_quantumL; % Preserve quantum numbers
H_hr.elementL = sc_elementL; % Preserve element information
H_hr.num = true; % Flag numerical storage
H_hr.coe = false; % Clear symbolic coefficients

end

% ==================== Helper Functions ====================

function H_hr = process_sparse(H_hr, OUT_WAN_NUM, repeatnum, fin_dir, glue_edges, vectorList)
%PROCESS_SPARSE Handle sparse matrix format Hamiltonian
    
    % Initialize sparse matrix array
    OUT_HnumL = repmat({sparse(OUT_WAN_NUM, OUT_WAN_NUM)}, 1, H_hr.NRPTS);
    WANNUM = H_hr.WAN_NUM;
    
    progress = CmdLineProgressBar('Processing sparse hoppings: ');
    for iN = 1:H_hr.NRPTS
        progress.print(iN, H_hr.NRPTS);
        ind_R = vectorList(iN,:);
        jump_fin = ind_R(fin_dir); % Finite direction hopping component
        
        % Extract non-zero elements
        [ilist, jlist, amplist] = find(H_hr.HnumL{iN});
        
        for ih = 1:length(amplist)
            i = ilist(ih);
            j = jlist(ih);
            amp = amplist(ih);
            
            % Apply repetitions along finite direction
            for icur_sc = 1:repeatnum
                hi = i + (icur_sc-1)*WANNUM;
                hj = j + (icur_sc + jump_fin - 1)*WANNUM;
                
                % Boundary condition handling
                [hj, valid] = apply_boundary_conditions(hj, OUT_WAN_NUM, glue_edges);
                if valid
                    OUT_HnumL{iN}(hi, hj) = amp;
                end
            end
        end
    end
    H_hr.HnumL = OUT_HnumL;
end

function H_hr = process_matrix(H_hr, OUT_WAN_NUM, repeatnum, fin_dir, glue_edges, vectorList, NRPT_seq)
%PROCESS_MATRIX Handle full matrix format Hamiltonian
    
    OUT_HnumL = zeros(OUT_WAN_NUM, OUT_WAN_NUM, H_hr.NRPTS);
    WANNUM = H_hr.WAN_NUM;
    NRPTS_record = H_hr.NRPTS;
    
    progress = CmdLineProgressBar('Processing matrix hoppings: ');
    for iN = 1:NRPTS_record
        progress.print(iN, H_hr.NRPTS);
        ind_R = vectorList(iN,:);
        jump_fin = ind_R(fin_dir);
        
        % Extract non-zero elements
        [ilist, jlist, amplist] = find(H_hr.HnumL(:,:,iN));
        
        for ih = 1:length(amplist)
            i = ilist(ih);
            j = jlist(ih);
            amp = amplist(ih);
            
            for icur_sc = 1:repeatnum
                hi = i + (icur_sc-1)*WANNUM;
                hj = j + (icur_sc + jump_fin - 1)*WANNUM;
                ind_R_tmp = ind_R;
                
                % Boundary condition handling
                [hj, valid, ind_R_tmp] = matrix_boundary_handling(...
                    hj, OUT_WAN_NUM, glue_edges, ind_R_tmp, fin_dir);
                
                if valid
                    % Update hopping vector index
                    [~, IH] = ismember(ind_R_tmp, vectorList, 'rows');
                    if IH == 0
                        IH = NRPTS_record + 1;
                        NRPTS_record = NRPTS_record + 1;
                        vectorList(IH,:) = ind_R_tmp;
                    end
                    NRPT_seq(IH) = true;
                    OUT_HnumL(hi, hj, IH) = amp;
                end
            end
        end
    end
    
    H_hr.vectorL = vectorList;
    H_hr.HnumL = OUT_HnumL(:,:,1:NRPTS_record);
    H_hr = H_hr.reseq(':', NRPT_seq); % Reindex hoppings
end

function H_hr = process_list(H_hr, OUT_WAN_NUM, repeatnum, fin_dir, glue_edges, vectorList)
%PROCESS_LIST Handle list format Hamiltonian
    
    OUT_HnumL = zeros(size(vectorList,1)*repeatnum, 1);
    OUT_vectorList = zeros(size(vectorList,1)*repeatnum, size(vectorList,2));
    WANNUM = H_hr.WAN_NUM;
    count = 0;
    
    progress = CmdLineProgressBar('Processing list hoppings: ');
    for ih = 1:H_hr.NRPTS
        progress.print(ih, H_hr.NRPTS);
        ind_R = vectorList(ih, 1:H_hr.Dim);
        jump_fin = ind_R(fin_dir);
        i = vectorList(ih, H_hr.Dim+1);
        j = vectorList(ih, H_hr.Dim+2);
        amp = H_hr.HnumL(ih);
        
        if abs(amp) > eps
            for icur_sc = 1:repeatnum
                hi = i + (icur_sc-1)*WANNUM;
                hj = j + (icur_sc + jump_fin - 1)*WANNUM;
                
                % Boundary condition handling
                [hj, valid] = apply_boundary_conditions(hj, OUT_WAN_NUM, glue_edges);
                
                if valid
                    count = count + 1;
                    OUT_HnumL(count) = amp;
                    OUT_vectorList(count,:) = [ind_R, hi, hj];
                end
            end
        end
    end
    
    % Trim unused preallocated space
    OUT_HnumL = OUT_HnumL(1:count);
    OUT_vectorList = OUT_vectorList(1:count,:);
    OUT_vectorList(:,fin_dir) = 0; % Flatten cutting direction
    
    H_hr.HnumL = OUT_HnumL;
    H_hr.vectorL = OUT_vectorList;
end

% ==================== Boundary Condition Utilities ====================

function [hj, valid] = apply_boundary_conditions(hj, OUT_WAN_NUM, glue_edges)
%APPLY_BOUNDARY_CONDITIONS Handle edge cases for hopping targets
    valid = true;
    if ~glue_edges
        if hj <= 0 || hj > OUT_WAN_NUM
            valid = false;
        end
    else
        hj = mod(hj-1, OUT_WAN_NUM) + 1; % 1-based indexing
    end
end

function [hj, valid, ind_R_tmp] = matrix_boundary_handling(hj, OUT_WAN_NUM, glue_edges, ind_R_tmp, fin_dir)
%MATRIX_BOUNDARY_HANDLING Special handling for matrix format with vector tracking
    valid = true;
    if ~glue_edges
        if hj <= 0 || hj > OUT_WAN_NUM
            valid = false;
        end
        ind_R_tmp(fin_dir) = 0;
    else
        ind_R_tmp(fin_dir) = floor((hj-1)/OUT_WAN_NUM);
        hj = mod(hj-1, OUT_WAN_NUM) + 1;
    end
end