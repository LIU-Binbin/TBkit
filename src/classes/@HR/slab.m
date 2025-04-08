function [EIGENCAR_slab,klist_l,kpoints_l,kpoints_name] = slab(H_hr,repeatnum,fin_dir,KPOINTS_slab,norb_enforce,fermi)
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
