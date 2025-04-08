function [nodes_s, nodes_r] = findnodes(Ham_obj, opts, kopts)
% Find the nodes (energy gap minima) in the Brillouin zone for a given Hamiltonian object.
% 
% Usage:
%   [nodes_s, nodes_r] = findnodes(Ham_obj, opts, kopts)
%
% Inputs:
%   Ham_obj      - Hamiltonian object (can be of types 'HR', 'Htrig', 'HK').
%   opts         - Options structure with fields:
%       .Num_Occupied  - Number of occupied bands (default: 0, i.e., half of bands).
%       .Gap_Threshold  - Energy gap threshold for detecting nodes (default: 0.0001).
%   kopts        - K-point options structure with fields:
%       .nk             - Number of k-points in each dimension [nk1, nk2, nk3].
%       .vk             - K-point grid vectors.
%       .original_point - Shift for the k-points to a different origin.
%       .mode           - Node search mode: 'corner' or 'center' (default: 'corner').
%
% Outputs:
%   nodes_s      - Nodes in fractional coordinates.
%   nodes_r      - Nodes in real-space coordinates.

arguments
    Ham_obj TBkit  % Hamiltonian object (must be of class 'TBkit').
    opts.Num_Occupied int8 = 0  % Default: 0 (half of the bands are occupied).
    opts.Gap_Threshold double = 0.0001  % Default threshold for gap.
    
    kopts.nk int8 = [10 10 10];  % Number of k-points in each direction [nk1, nk2, nk3].
    kopts.vk = [1 0 0; 0 1 0; 0 0 1];  % K-point grid vectors (default identity matrix).
    kopts.original_point = [-0.5 -0.5 -0.5];  % Original point to shift the k-points.
    kopts.mode {mustBeMember(kopts.mode, {'corner', 'center'})} = 'corner';  % Node search mode.
    kopts.edge {mustBeMember(kopts.edge, {'half', 'full'})} = "half";  % Edge specification for nodes.
end

% Convert named arguments to cell for use in functions
koptscell = namedargs2cell(kopts);

% Set number of occupied bands
if opts.Num_Occupied == 0
    nbands = Ham_obj.Nbands;  % Total number of bands
    occu = nbands / 2;        % Default: half of the bands are occupied
else
    occu = opts.Num_Occupied;
end

% Generate k-point meshes (both fractional and real-space)
[klist_s, klist_r] = kmesh_gen(Ham_obj, [], koptscell{:});

% Define function to calculate the energy gap at a given k-point
switch class(Ham_obj)
    case 'HR'
        get_gap_kpoint = @(kpoint) get_gap(Ham_obj, kpoint, occu);
    case {'Htrig', 'HK'}
        Hfun = Ham_obj.Hfun;  % Handle for Hamiltonian function
        get_gap_kpoint = @(kpoint) get_gap_Hfun(Hfun, kpoint, occu);
end

nkpts = size(klist_s, 1);  % Total number of k-points
nodes_tmp = [];            % Temporary array to store nodes

% Search for nodes with energy gap below threshold using optimization
pb = TBkit_tool_outer.CmdLineProgressBar('FindNode: ');  % Progress bar for visualization
count = 0;  % Counter for the number of nodes found

% Loop over all k-points to find nodes
for i = 1:nkpts
    % Minimize the energy gap at each k-point using fminsearch
    [kout, gapout] = fminsearch(get_gap_kpoint, klist_s(i, :));
    
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
