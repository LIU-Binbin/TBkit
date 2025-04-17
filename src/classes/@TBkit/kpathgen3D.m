function TBkitobj = kpathgen3D(TBkitobj, KPOINTS_name, nodes, kpoints_name_tmp)
%KPATHGEN3D Generate k-point paths for band structure calculations and plotting
%   TBkitobj = KPATHGEN3D(TBkitobj, KPOINTS_name, nodes, kpoints_name_tmp) 
%   creates k-point paths for quantum mechanical calculations using either
%   KPOINTS file input or direct coordinate input. The function handles both
%   fractional and Cartesian coordinates, and stores results in TBkitobj.
%
%   Inputs:
%       TBkitobj        - Tight-binding toolkit object containing lattice info
%       KPOINTS_name    - (Optional) Name of KPOINTS file or k-points matrix. 
%                         Default: 'KPOINTS'
%       nodes           - (Optional) Number of nodes per path segment. If empty, 
%                         reads from KPOINTS file
%       kpoints_name_tmp- (Optional) Names of high-symmetry points
%
%   Output:
%       TBkitobj - Updated object containing:
%           klist_cart   : Cartesian coordinates of k-points (Nx3 matrix)
%           klist_frac   : Fractional coordinates of k-points (Nx3 matrix) 
%           klist_l      : Normalized path length for plotting (Nx1 vector)
%           kpoints_l    : High-symmetry point positions along path (Mx1)
%           kpoints_frac : Fractional coordinates of high-symmetry points (Mx3)
%           kpoints_name : Names of high-symmetry points (cell array)
%
%   Usage Examples:
%       1. Default KPOINTS file usage:
%          TBkitobj = kpathgen3D(TBkitobj);
%
%       2. Custom k-points matrix with 50 nodes per segment:
%          kpoints = [0 0 0; 0.5 0 0; 0.5 0.5 0];
%          TBkitobj = kpathgen3D(TBkitobj, kpoints, 50, {'Γ','X','M'});
%
%    Supported File Formats:
%       KPOINTS file must use 'Line-Mode' and 'Reciprocal' coordinates
%
%    Change Log:
%       2020-12-03 - Initial implementation
%       2022-11-06 - Updated argument handling
%
%    Author: parkman <parkman@buaa.edu.cn>

arguments
    TBkitobj
    KPOINTS_name = 'KPOINTS';
    nodes = [];
    kpoints_name_tmp = [];
end

% Initialize lattice vectors if empty
if isempty(TBkitobj.Rm)
    TBkitobj = TBkitobj.input_Rm();
end

% Handle different input modes
if isempty(nodes)
    % Read from KPOINTS file
    [kpoints, nodes, kpoints_name_tmp] = KPOINTS_read(KPOINTS_name);
else
    % Use direct k-points matrix input
    kpoints = KPOINTS_name;
end

% Generate k-path using reciprocal lattice vectors
[TBkitobj.klist_cart, TBkitobj.klist_frac, TBkitobj.klist_l, ...
 TBkitobj.kpoints_l, TBkitobj.kpoints_frac] = ...
    TBkit.kpathgen(kpoints, nodes, TBkitobj.Gk);

% Store high-symmetry point labels
TBkitobj.kpoints_name = kpoints_name_tmp;
end