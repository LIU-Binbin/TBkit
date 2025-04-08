function [DOSCAR_l,DOSCAR_b,DOSCAR_r,w_list,klist_l,kpoints_l,kpoints_name] = surf(H_hr,w_range,fin_dir,KPOINTS_surf,principle_layer,eta,fermi,mode)
if nargin < 8
mode = 'Green_iter';
end
if nargin < 7
fermi = 0;
end
if nargin < 6
eta = 0.01;
end
if nargin < 5
principle_layer = 2;
end
if nargin < 4
KPOINTS_surf = 'KPOINTS_surf';
end
if nargin < 3
fin_dir = 2 ;
end
if nargin < 2
w_range = [-1,1,100];
end
w_list = linspace(w_range(1),w_range(2),w_range(3))+fermi;
H_hr_surf = H_hr < KPOINTS_surf;
[H00_H01_cell_list_1,H00_H01_cell_list_2] = H_hr_surf.H00_H11_cell_list_gen(fin_dir,principle_layer);
GREENCAR = HR.GREENCAR_gen(w_list,eta,H00_H01_cell_list_1,H00_H01_cell_list_2,mode);
DOSCAR_b = HR.DOSCAR_gen(GREENCAR.bulk,'green');
DOSCAR_l = HR.DOSCAR_gen(GREENCAR.surf_l,'green');
DOSCAR_r = HR.DOSCAR_gen(GREENCAR.surf_r,'green');
[klist_l,kpoints_l,kpoints_name] = H_hr_surf.kpath_information();
w_list = w_list - fermi;
end
