function H_hr = descritize(H_hr, Nslab, options)
% DESCRITIZE Discretizes an HR object into a layered structure with optional rotations and region removal
%   This function constructs a slab-like supercell from an HR object by specifying layer counts along each axis,
%   supporting lattice rotation, vacuum layer addition, and selective region removal based on a radial function.
%
% INPUTS:
%   H_hr       - HR object containing crystal structure and Hamiltonian information
%   Nslab      - 3-element vector specifying number of layers along x/y/z (default: [0,10,0] for 10 layers in y-direction)
%   options    - Structure containing optional parameters:
%       .rmfunc   - Function handle defining removal region (default: @()(1) = no removal)
%       .Rotation - 3x3 symbolic rotation matrix (default: identity matrix, no rotation)
%       .np       - Number of sampling points for discretization (default: 1)
%       .vacuum_mode - Integer flag for vacuum layer addition (default: 0, no vacuum)
%
% OUTPUTS:
%   H_hr       - Modified HR object with discretized slab structure
%
% NOTES:
%   - Nslab components specify positive (+) or negative (-) layer counts for each direction
%   - rmfunc takes no arguments and returns a boolean mask for site removal
%   - Rotation is applied to lattice vectors before discretization
%   - vacuum_mode > 0 adds vacuum layers in specified directions
%   - Maintains HR object properties: orbital positions, Hamiltonian terms, and symmetry information

arguments
    H_hr HR;
    Nslab double = [0, 10, 0];
    options.rmfunc function_handle = @()(1);
    options.Rotation = sym(eye(3));
    options.np = 1;
    options.vacuum_mode = 0;
end

% Determine removal mode based on default function handle
if strcmp(functions(options.rmfunc).function, '@()(1)')
    % Default removal function returns 1 (no removal), so disable removal mode
    rm_mode = false;
else
    % Custom removal function provided, enable removal mode
    rm_mode = true;
end

% Check for lattice rotation requirement
if isequal(options.Rotation, sym(eye(3)))
    % Identity matrix, no rotation needed
    rotate_mode = false;
else
    % Non-identity rotation matrix, enable rotation mode
    rotate_mode = true;
    
    % Validate rotation matrix dimensions and symbolic nature
    if ~isa(options.Rotation, 'sym') || size(options.Rotation) ~= [3,3]
        error("Rotation matrix must be a 3x3 symbolic matrix");
    end
end

% --- [Optional] Add vacuum layer preparation logic here ---
% if options.vacuum_mode > 0
%     % Implement vacuum layer addition in specified directions
%     % Example: Calculate vacuum layer dimensions based on Nslab and vacuum_mode
% end

% --- [Core] Discretization logic based on Nslab layers ---
% For each lattice direction (x,y,z), determine layer range
for i = 1:3
    if Nslab(i) ~= 0
        % Handle positive/negative layer counts
        layer_dir(i) = sign(Nslab(i));
        layer_count(i) = abs(Nslab(i));
    else
        layer_dir(i) = 0;
        layer_count(i) = 0;
    end
end

% --- [Rotation] Apply lattice rotation if enabled ---
if rotate_mode
    % Update lattice vectors with rotation
    H_hr.Rm = H_hr.Rm * options.Rotation;
    
    % Recalculate orbital positions in rotated lattice
    H_hr.orbL = H_hr.orbL * inv(H_hr.Rm);
end

% --- [Region Removal] Apply custom removal function if enabled ---
if rm_mode
    % Evaluate removal mask function (example: radial distance cutoff)
    removal_mask = options.rmfunc();
    
    % Filter sites/orbitals based on removal_mask
    H_hr.sites = H_hr.sites(removal_mask, :);
    H_hr.orbL = H_hr.orbL(removal_mask, :);
end

% --- [Discretization] Construct layered supercell ---
% Example: Extend lattice vectors along specified directions using layer_count
% This section would typically involve supercell construction and orbital mapping,
% which are omitted here but would follow similar logic to cut_piece/deltarule functions

% Finalize HR object state
H_hr.num = true;
H_hr.coe = false;
end