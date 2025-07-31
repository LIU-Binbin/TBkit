function [nodes_s, nodes_r] = findnodes(Ham_obj, options)
arguments
    Ham_obj TBkit  % Hamiltonian object (must be of class 'TBkit').
    options.Num_Occupied int8 = 0  % Default: 0 (half of the bands are occupied).
    options.Gap_Threshold double = 0.0001  % Default threshold for gap.
    
    options.kstart(1,3) double = [0 0 0]
    options.kdir1 (1,3) double = [1 0 0]
    options.kdir2 (1,3) double = [0 1 0]
    options.kdir3 (1,3) double = [0 0 1]
    options.Nk1 double = 10
    options.Nk2 double = 10
    options.Nk3 double = 10
end

% Convert named arguments to cell for use in functions
optscell = namedargs2cell(options);

% Set number of occupied bands
if opts.Num_Occupied == 0
    nbands = Ham_obj.Nbands;  % Total number of bands
    occu = nbands / 2;        % Default: half of the bands are occupied
else
    occu = opts.Num_Occupied;
end

[~, klist_frac] = kmeshgen(Ham_obj.Rm, optscell{:});

% Define function to calculate the energy gap at a given k-point
switch class(Ham_obj)
    case 'HR'
        get_gap_kpoint = @(kpoint) get_gap(Ham_obj, kpoint, occu);
    case {'Htrig', 'HK'}
        Hfun = Ham_obj.Hfun;  % Handle for Hamiltonian function
        get_gap_kpoint = @(kpoint) get_gap_Hfun(Hfun, kpoint, occu);
end

nkpts = size(klist_frac, 1);  % Total number of k-points
nodes_tmp = [];            % Temporary array to store nodes

% Search for nodes with energy gap below threshold using optimization
pb = CmdLineProgressBar('FindNode: ');  % Progress bar for visualization
count = 0;  % Counter for the number of nodes found

% Loop over all k-points to find nodes
for i = 1:nkpts
    % Minimize the energy gap at each k-point using fminsearch
    [kout, gapout] = fminsearch(get_gap_kpoint, klist_frac(i, :));
    
    % If energy gap is smaller than the threshold, consider it a node
    if gapout < opts.Gap_Threshold
        nodes_tmp = [nodes_tmp; kout];  % Add node to temporary list
        count = count + 1;
    end
    
    % Update progress bar
    msgtail = [' hit:', num2str(count), '/', num2str(nkpts)];
    pb.print(i, nkpts, msgtail);
end
pb.delete();  % Delete progress bar when done

% If no nodes found, return empty arrays
if isempty(nodes_tmp)
    nodes_s = [];
    nodes_r = [];
    return;
end

% Ensure k-points are within the valid Brillouin zone (if nk = 1 for any direction)
for i = 1:3
    if kopts.nk(i) == 1
        nodes_tmp(:, i) = 0;  % Set k-point components to 0 if nk = 1
    end
end

% Convert k-points to fractional coordinates
switch class(Ham_obj)
    case "HR"
        nodes_s = nodes_tmp;  % Use nodes directly for 'HR' type Hamiltonians
    case {"Htrig", "HK"}
        nodes_s = nodes_tmp / Ham_obj.Gk;  % Normalize by the reciprocal lattice vectors
end

% Shift k-points to a user-defined block
nodes_s = kshift(nodes_s, [kopts.original_point; kopts.vk]);  % Apply the user-defined shift
nodes_s = uniquetol(nodes_s, 1e-4, 'ByRows', true);  % Remove duplicates with tolerance

% Convert fractional k-points back to real-space coordinates
nodes_r = nodes_s * Ham_obj.Gk;  % Multiply by reciprocal lattice vectors to get real-space coordinates

end


% Helper function to calculate the energy gap for 'Htrig' or 'HK' type Hamiltonians
function gap = get_gap_Hfun(Hfun, kpoint, occu)
    H = Hfun(kpoint(1), kpoint(2), kpoint(3));  % Get Hamiltonian matrix at k-point
    EIGENCAR = eig((H + H') / 2);  % Diagonalize the Hamiltonian
    gap = EIGENCAR(occu + 1, :) - EIGENCAR(occu, :);  % Calculate the gap between bands
end

% Helper function to calculate the energy gap for 'HR' type Hamiltonians
function gap = get_gap(Ham_obj, kpoint, occu)
    EIGENCAR = Ham_obj.EIGENCAR_gen('klist', kpoint, 'printmode', false);  % Generate eigenvalues at k-point
    gap = EIGENCAR(occu + 1, :) - EIGENCAR(occu, :);  % Calculate the gap between bands
end
