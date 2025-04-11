function [DOSCAR_l,DOSCAR_b,DOSCAR_r,w_list,klist_l,kpoints_l,kpoints_name] = surf(H_hr,w_range,fin_dir,KPOINTS_surf,principle_layer,eta,fermi,mode)
% SURF Calculate surface spectral properties
%
%   [DOSCAR_L,DOSCAR_B,DOSCAR_R,WLIST,KLIST_L,KPOINTS_L,KPOINTS_NAME] = SURF(...)
%   computes surface density of states and spectral properties
%
%   Inputs:
%       H_hr - Bulk HR object
%       w_range - Energy range [default: [-1,1,100]]
%       fin_dir - Fin direction [default: 2]
%       KPOINTS_surf - K-points file [default: 'KPOINTS_surf']
%       principle_layer - Principle layers [default: 2]
%       eta - Broadening [default: 0.01]
%       fermi - Fermi level [default: 0]
%       mode - Calculation mode [default: 'Green_iter']
%   Outputs:
%       DOSCAR_l - Left surface DOS
%       DOSCAR_b - Bulk DOS
%       DOSCAR_r - Right surface DOS
%       w_list - Energy list
%       klist_l - K-point list
%       kpoints_l - K-point coordinates
%       kpoints_name - K-point names
%
%   Notes:
%       - Uses Green's function methods
%       - Handles both left and right surfaces
%       - Returns k-path information
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
