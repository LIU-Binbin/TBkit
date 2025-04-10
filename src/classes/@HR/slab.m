function [EIGENCAR_slab,klist_l,kpoints_l,kpoints_name] = slab(H_hr,repeatnum,fin_dir,KPOINTS_slab,norb_enforce,fermi)
% SLAB Construct slab system and calculate band structure
%
%   [EIGENCAR,KLIST_L,KPOINTS_L,KPOINTS_NAME] = SLAB(H_HR,REPEATNUM,FIN_DIR,KPOINTS_SLAB,NORB_ENFORCE,FERMI)
%   creates slab system and calculates band structure
%
%   Inputs:
%       H_hr - Bulk HR object
%       repeatnum - Slab repetition number [default: 10]
%       fin_dir - Fin direction [default: 2]
%       KPOINTS_slab - K-points file [default: 'KPOINTS_slab']
%       norb_enforce - Orbital enforcement [default: -1]
%       fermi - Fermi level [default: 0]
%   Outputs:
%       EIGENCAR_slab - Eigenvalue array
%       klist_l - K-point list
%       kpoints_l - K-point coordinates
%       kpoints_name - K-point names
%
%   Notes:
%       - Uses cut_piece to create slab
%       - Calculates band structure along specified path
%       - Returns eigenvalues and k-space information
if nargin < 6
fermi = 0;
end
if nargin < 5
norb_enforce  = -1;
end
if nargin <4
KPOINTS_slab = 'KPOINTS_slab';
end
if nargin < 3
fin_dir     =  2;
end
if nargin < 2
repeatnum   = 10;
end
glue_edges  = false;
vacuum_mode = 1;
H_hr_slab = H_hr.cut_piece(repeatnum,fin_dir,glue_edges,vacuum_mode);
H_hr_slab = H_hr_slab < KPOINTS_slab;
EIGENCAR_slab = H_hr_slab.EIGENCAR_gen('fermi',fermi,'norb',norb_enforce);
[klist_l,kpoints_l,kpoints_name] = H_hr_slab.kpath_information();
end
