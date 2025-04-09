function H_hr_n = Connect2D(H_hr_n, H_hr_unitcell, opt)
% CONNECT2D Construct 2D lattice from unit cell Hamiltonians with various connection schemes
%
%   H_hr_n = CONNECT2D(H_hr_n, H_hr_unitcell, opt) builds a 2D extended
%   Hamiltonian by connecting multiple copies of a unit cell Hamiltonian
%   according to specified connection patterns and symmetries.
%
%   Inputs:
%       H_hr_n         - Target HR object for extended system
%       H_hr_unitcell  - Unit cell HR object to replicate
%       opt            - Connection options structure:
%           Sequ          : Unit cell replication sequence [1:20]
%           findir        : Principal connection direction (1=x, 2=y)
%           Subsequ       : Connection mode ('normal','reverse','combined')
%           reverse_sequ  : Orbital reversal sequence [default: reverse WAN_NUM]
%           combined_norm : Normal sequence indices for combined mode [2:19]
%
%   Output:
%       H_hr_n - Extended HR object with connected 2D system
%
%   Features:
%       - Supports normal, reversed, and mixed connection schemes
%       - Maintains Hermiticity through conjugate hopping terms
%       - Handles 2D connection patterns in both x and y directions
%       - Includes symmetry checks for reversal operations
%
%   Example:
%       opt = struct('Sequ',1:5,'findir',2,'Subsequ','combined');
%       H_hr_2D = Connect2D(H_hr_empty, H_hr_unit, opt);
%
%   See also HR, SET_HOP, AUTOHERMI

    arguments
        H_hr_n HR
        H_hr_unitcell HR
        opt.Sequ (1,:) double {mustBeInteger,mustBePositive} = 1:20
        opt.findir (1,1) double {mustBeMember(opt.findir,[1,2])} = 2
        opt.Subsequ {mustBeMember(opt.Subsequ,{'normal','reverse','combined'})} = 'normal'
        opt.reverse_sequ (1,:) double = H_hr_unitcell.WAN_NUM:-1:1
        opt.combined_norm (1,:) double {mustBeInteger,mustBePositive} = 2:19
    end

    % Validate system dimensions
    if H_hr_unitcell.Dim ~= 2 || H_hr_n.Dim ~= 2
        error('Connect2D requires 2D HR objects');
    end
    
    % Initialize system parameters
    WANNUM = H_hr_n.WAN_NUM;
    WANNUM_unit = H_hr_unitcell.WAN_NUM;
    
    % Calculate replication parameters
    if opt.findir == 2
        process_y_direction();
    else
        process_x_direction();
    end

    % Nested processing functions
    function process_y_direction()
        Nrep_y = length(opt.Sequ);
        Nrep_x = WANNUM/WANNUM_unit/Nrep_y;
        validate_replication(Nrep_x);
        
        % Process nearest neighbor connections
        for k = -1:1
            vec_pattern = create_vector_pattern(1, k); % [1, k] for y-direction
            process_connections(vec_pattern, [1 0], Nrep_x, Nrep_y, k);
        end
    end

    function process_x_direction()
        Nrep_x = length(opt.Sequ);
        Nrep_y = WANNUM/WANNUM_unit/Nrep_x;
        validate_replication(Nrep_y);
        
        % Process nearest neighbor connections
        for k = -1:1
            vec_pattern = create_vector_pattern(k, 1); % [k, 1] for x-direction
            process_connections(vec_pattern, [0 1], Nrep_x, Nrep_y, k);
        end
    end

    % Helper functions
    function validate_replication(rep_factor)
        if mod(rep_factor,1) ~= 0
            error('Invalid WAN_NUM: Must be integer multiple of unit cells');
        end
    end

    function vec = create_vector_pattern(vx, vy)
        % Create 2D vector pattern based on connection direction
        vec = [vx, vy];
    end

    function process_connections(vec_pattern, base_vector, Nrep_main, Nrep_ortho, k)
        % Main connection processing engine
        slct = ismember(H_hr_unitcell.vectorL(:,1:2), vec_pattern, 'rows');
        if ~any(slct), return; end
        
        % Extract hopping parameters
        [row, col, num, num_rev] = get_hopping_parameters(slct);
        
        % Validate reversal symmetry if needed
        if contains(opt.Subsequ, 'reverse') && ~check_reversal_symmetry(num, num_rev)
            warning('Potential symmetry violation in unit cell hoppings');
        end
        
        % Apply connection pattern
        for i = 1:length(opt.Sequ)
            [labR, labL] = calculate_cell_indices(Nrep_main, Nrep_ortho, i, k);
            if labL <= 0 || labL > Nrep_ortho, continue; end
            
            apply_hopping_terms(row, col, num, num_rev, labR, labL, base_vector);
        end
    end

    function [row, col, num, num_rev] = get_hopping_parameters(slct)
        % Extract hopping parameters based on subsequence mode
        HnumL_tmp = H_hr_unitcell.HnumL(:,:,slct);
        
        switch opt.Subsequ
            case 'normal'
                [row, col, num] = find(HnumL_tmp);
                num_rev = [];
                
            case 'reverse'
                HnumL_reverse = HnumL_tmp(:, opt.reverse_sequ);
                [row, col, num] = find(HnumL_reverse);
                num_rev = flip(num);
                
            case 'combined'
                [row_norm, col_norm, num_norm] = find(HnumL_tmp);
                HnumL_reverse = HnumL_tmp(:, opt.reverse_sequ);
                [row_rev, col_rev, num_rev] = find(HnumL_reverse);
                num = num_norm;
                row = [row_norm; row_rev];
                col = [col_norm; col_rev];
        end
    end

    function valid = check_reversal_symmetry(num, num_rev)
        % Validate reflection symmetry for reversal operations
        if isempty(num_rev), valid = true; return; end
        valid = all(abs(num - num_rev) < 1e-6);
    end

    function [labR, labL] = calculate_cell_indices(Nrep_main, Nrep_ortho, i, k)
        % Calculate cell indices based on connection direction
        if opt.findir == 2
            labR = Nrep_ortho*(Nrep_main-1) + i;
            labL = opt.Sequ(i) + k;
        else
            labR = i*Nrep_ortho;
            tmp = (opt.Sequ(i)-1+k);
            labL = Nrep_ortho*tmp + 1;
        end
    end

    function apply_hopping_terms(row, col, num, num_rev, labR, labL, vector)
        % Apply hopping terms with Hermitian conjugation
        if strcmp(opt.Subsequ, 'combined')
            apply_combined_hopping(labR, labL, vector);
        else
            for j = 1:length(row)
                idxR = WANNUM_unit*(labR-1) + row(j);
                idxL = WANNUM_unit*(labL-1) + col(j);
                
                % Add forward and reverse hoppings
                H_hr_n = H_hr_n.set_hop(num(j), idxR, idxL, vector, 'add');
                H_hr_n = H_hr_n.set_hop(conj(num(j)), idxL, idxR, -vector, 'add');
            end
        end
    end

    function apply_combined_hopping(labR, labL, vector)
        % Special handling for combined normal/reverse connections
        persistent count; if isempty(count), count = 1; end
        
        if count <= length(opt.combined_norm) && ismember(labR, opt.combined_norm)
            % Apply normal hopping terms
            for j = 1:length(row_norm)
                idxR = WANNUM_unit*(labR-1) + row_norm(j);
                idxL = WANNUM_unit*(labL-1) + col_norm(j);
                H_hr_n = H_hr_n.set_hop(num_norm(j), idxR, idxL, vector, 'add');
                H_hr_n = H_hr_n.set_hop(conj(num_norm(j)), idxL, idxR, -vector, 'add');
            end
            count = count + 1;
        else
            % Apply reversed hopping terms
            for j = 1:length(row_rev)
                idxR = WANNUM_unit*(labR-1) + row_rev(j);
                idxL = WANNUM_unit*(labL-1) + col_rev(j);
                H_hr_n = H_hr_n.set_hop(num_rev(j), idxR, idxL, vector, 'add');
                H_hr_n = H_hr_n.set_hop(conj(num_rev(j)), idxL, idxR, -vector, 'add');
            end
        end
    end
end