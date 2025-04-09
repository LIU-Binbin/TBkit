function H_hr = cut_orb(H_hr, rm_list, options)
% CUT_ORB Remove specified orbitals from HR object
%
%   H_hr = CUT_ORB(H_hr, rm_list, options) removes orbitals from an HR object
%   based on either direct indices or a selection function. Maintains
%   Hamiltonian consistency after orbital removal.
%
%   Inputs:
%       H_hr     - HR object containing orbital information
%       rm_list  - Orbital removal specification:
%                  - Vector of indices to remove
%                  - Logical mask (true = remove)
%                  - Empty (use removal function)
%       options  - Optional parameters:
%           rmfunc - Function handle for dynamic removal (@()(1))
%
%   Output:
%       H_hr - Modified HR object with specified orbitals removed
%
%   Features:
%       - Automatic index validation
%       - Supports both direct and functional orbital selection
%       - Maintains Hamiltonian consistency
%
%   Example:
%       % Remove orbitals 3 and 5 directly
%       H_hr = cut_orb(H_hr, [3,5]);
%
%       % Remove orbitals using custom function (e.g., z > 0)
%       opt.rmfunc = @(x,y,z) z > 0;
%       H_hr = cut_orb(H_hr, [], opt);
%
%   See also HR, RESEQ, SELECT_ORB

    arguments
        H_hr HR
        rm_list {mustBeVector} = []
        options.rmfunc function_handle = @(x,y,z) false(size(x))
    end

    % Extract orbital coordinates
    orb_tmp = H_hr.orbL;
    num_orbitals = size(orb_tmp, 1);
    
    % Generate removal mask
    if isempty(rm_list)
        % Use removal function if provided
        if ~isequal(options.rmfunc, @(x,y,z) false(size(x)))
            validate_rmfunc(options.rmfunc, num_orbitals);
            rm_mask = options.rmfunc(orb_tmp(:,1), orb_tmp(:,2), orb_tmp(:,3));
        else
            % Default: keep all orbitals
            rm_mask = false(num_orbitals, 1);
        end
    else
        % Convert input to logical mask
        rm_mask = process_rmlist(rm_list, num_orbitals);
    end
    
    % Validate final removal mask
    if all(rm_mask)
        error('Cannot remove all orbitals from HR object');
    end
    
    % Create keep mask and resequence
    keep_mask = ~rm_mask;
    H_hr = H_hr.reseq(keep_mask);

    % Nested helper functions
    function rm_mask = process_rmlist(input_list, max_index)
        % Convert various input types to logical mask
        if islogical(input_list)
            validate_logical_mask(input_list, max_index);
            rm_mask = input_list;
        else
            validate_indices(input_list, max_index);
            rm_mask = false(max_index, 1);
            rm_mask(input_list) = true;
        end
    end

    function validate_rmfunc(func_handle, expected_size)
        % Test rmfunc output dimensions
        test_output = func_handle(0,0,0); % Single coordinate test
        if ~islogical(test_output) || numel(test_output) ~= 1
            error('rmfunc must return logical values matching orbital count');
        end
        
        % Full validation
        full_output = func_handle(orb_tmp(:,1), orb_tmp(:,2), orb_tmp(:,3));
        if ~islogical(full_output) || numel(full_output) ~= num_orbitals
            error('rmfunc must return logical vector matching orbital count');
        end
    end

    function validate_logical_mask(mask, expected_length)
        if numel(mask) ~= expected_length
            error('Logical mask length (%d) must match orbital count (%d)',...
                  numel(mask), expected_length);
        end
    end

    function validate_indices(indices, max_index)
        if any(indices < 1) || any(indices > max_index)
            error('Orbital indices out of range [1-%d]', max_index);
        end
        if numel(unique(indices)) ~= numel(indices)
            error('Duplicate indices in removal list');
        end
    end
end