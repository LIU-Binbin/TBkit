function TBkitobj = nn(TBkitobj, search_range, Accuracy, Rlength_cut, options)
%% Calculate nearest neighbor information for primitive cell structure
% This function computes neighbor relationships between atomic orbitals in
% a crystal structure based on lattice vectors and orbital positions.
%
% Inputs:
%   TBkitobj   - TBkit object containing structural information:
%                  Rm (3x3 lattice vectors), orbL (Nx3 orbital positions)
%   search_range - [Optional] 3-element vector specifying search range
%                  in lattice units (default: [0 0 0])
%   Accuracy     - [Optional] Numerical tolerance for distance comparison
%                  (default: 1e-4)
%   Rlength_cut  - [Optional] Maximum neighbor distance cutoff in Angstroms
%                  (default: 5)
%   options      - [Optional] Structure with additional parameters:
%       MAX_NN_LENGTH - Maximum allowed neighbor pairs (default: 1e7)
%       onsite        - Include onsite terms flag (default: false)
%
% Output:
%   TBkitobj   - Updated object with neighbor information:
%                  nn_store (neighbor list), Rnn (unique displacement vectors)
%
% Features:
% - Supports symbolic computation mode when Rm contains symbolic expressions
% - Progress tracking for large systems
% - Space-efficient sparse storage format
% - Automatic unique vector identification and sorting

%% Argument validation and initialization
arguments
    TBkitobj
    search_range = zeros(1, TBkitobj.Dim)
    Accuracy = 1e-4
    Rlength_cut = 5
    options.MAX_NN_LENGTH = 10000000
    options.onsite = false
end

%% Validate orbital data existence
if isempty(TBkitobj.orbL)
    error('Orbital positions (orbL) not loaded. Load structure data first.');
end

%% Handle precision warnings
if Accuracy > 1
    warning('High accuracy threshold detected. Verify units/scale.');
end

%% Initialize system parameters
sites_num = size(TBkitobj.orbL, 1);
[search_rangex, search_rangey, search_rangez] = deal(search_range(1), search_range(2), search_range(3));
sym_mode = isa(TBkitobj.Rm, 'sym');

%% Configure computational mode
if sym_mode
    Accuracy = -1;  % Disable numerical checks for symbolic mode
    Rm_ = TBkitobj.Rm;
else
    Rm_ = double(TBkitobj.Rm);
end

%% Initialize sparse storage
max_nn_length = min(...
    sites_num^2 * prod(search_range+1)^2, ...
    options.MAX_NN_LENGTH);
nn_sparse = zeros(max_nn_length, 10);  % Preallocate neighbor matrix
Rnn_list = [];                         % Initialize displacement list
count  = 1;
%% Main neighbor calculation loop
progress = CmdLineProgressBar('Neighbor Search Progress: ');
for j = 1:sites_num
    orb_j = double(TBkitobj.orbL(j,:));

    for i = 1:sites_num
        orb_i = double(TBkitobj.orbL(i,:));

        % Generate neighbor pairs for current orbital combination
        [nn_sparse_temp, displacements] = TBkitobj.nn_sparse_gen(...
            orb_i, orb_j, Rm_, search_rangex, search_rangey, search_rangez,...
            Accuracy, Rlength_cut, options.onsite);

        % Store results with orbital indices
        if ~isempty(nn_sparse_temp)
            countadd = size(nn_sparse_temp,1);
            nn_sparse_temp(:,1) = i;
            nn_sparse_temp(:,2) = j;
            nn_sparse(count:count+countadd-1,:) = nn_sparse_temp;
            count = count+countadd;
            Rnn_list = [Rnn_list; displacements];
        end
    end
    progress.print(j, sites_num, ' Processing orbitals...');
end
progress.delete();

%% Post-processing and cleanup
nn_sparse(count:max_nn_length,:) = [];% Trim unused preallocated space

%% Generate unique displacement vectors
TBkitobj.Rnn = sort(unique(Rnn_list, 'rows', 'stable'));

%% Assign neighbor levels with optional onsite exclusion
level_offset = double(logical(options.onsite));  % Adjust level numbering
progress = CmdLineProgressBar('Level Assignment: ');
for idx = 1:size(nn_sparse,1)
    [~, level] = ismember(nn_sparse(idx,9), TBkitobj.Rnn);
    nn_sparse(idx,10) = level - level_offset;
    progress.print(idx, size(nn_sparse,1), ' Processing bonds...');
end
progress.delete();

%% Finalize output storage
TBkitobj.nn_store = nn_sparse;
end