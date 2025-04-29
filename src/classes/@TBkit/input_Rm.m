function TBkitobj = input_Rm(TBkitobj,Rm)
%INPUT_RM Set lattice vectors for TBkit object
%   TBkitobj = input_Rm(TBkitobj, Rm) assigns lattice vectors to a TBkit object
%
%   Inputs:
%       TBkitobj - Target TBkit object
%       Rm       - Lattice vectors (3x3 matrix) or POSCAR filename
%
%   Output:
%       TBkitobj - Updated object with Rm property set
%
%   If no Rm provided, attempts to read from POSCAR file
%
%   Example:
%       TB = input_Rm(TB, [1 0 0; 0 1 0; 0 0 1]);
if nargin <2 && exist('POSCAR','file')
    [Rm,~,~,~,~] = POSCAR_read('POSCAR','vasp');
elseif nargin <2 && ~exist('POSCAR','file')
    Rm = eye(TBkitobj.Dim);
    %warning('POSCAR or Rm needed');
else
    if isa(Rm,'string') || isa(Rm,'char')
        [Rm,~,~,~,~] = POSCAR_read(Rm,'vasp');
    end
end
Rm = (Rm);
TBkitobj.Rm = Rm;
end