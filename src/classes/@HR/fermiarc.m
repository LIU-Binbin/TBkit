function [DOSCAR_l,DOSCAR_b,DOSCAR_r,klist1,klist2] = fermiarc(H_hr,w_arc,fin_dir,kmesh,kfermiarc,principle_layer,eta,fermi,mode)
% FERMIARC Calculate Fermi arc surface states
%
%   [DOSCAR_l,DOSCAR_b,DOSCAR_r,klist1,klist2] = FERMIARC(H_hr,w_arc,fin_dir,kmesh,kfermiarc,principle_layer,eta,fermi,mode)
%   computes surface Green's functions and density of states for Fermi arcs.
%
%   INPUT ARGUMENTS:
%       H_hr - Bulk Hamiltonian in HR format
%       w_arc - Energy range (relative to fermi)
%       fin_dir - Termination direction (1-3, default: 3)
%       kmesh - k-point mesh dimensions (default: [100 100])
%       kfermiarc - k-path basis (default: standard square path)
%       principle_layer - Number of principle layers (default: 2)
%       eta - Broadening parameter (default: 0.01)
%       fermi - Fermi level (default: 0)
%       mode - Green's function calculation mode (default: 'Green_iter')
%
%   OUTPUT ARGUMENTS:
%       DOSCAR_l - Left surface DOS
%       DOSCAR_b - Bulk DOS
%       DOSCAR_r - Right surface DOS
%       klist1,klist2 - k-point coordinates
%
%   NOTES:
%       - Uses recursive Green's function method
%       - Returns reshaped DOSCARs matching kmesh dimensions
%
%   SEE ALSO:
%       HR, GREENCAR_gen, DOSCAR_gen
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

if nargin < 9
mode = 'Green_iter';
end
if nargin < 8
fermi = 0;
end
if nargin < 7
eta = 0.01;
end
if nargin < 6
principle_layer = 2;
end
if nargin < 5
kfermiarc = [-0.5, -0.5, 0;...
1.0 , 0,   0;...
0   , 1.0,   0;...
];
end
if nargin < 4
kmesh = [100 100];
end
if nargin < 3
fin_dir = 3 ;
end
if nargin < 2
w_arc = 0;
end
w_arc = w_arc + fermi;
H_hr_arc = H_hr;
[H_hr_arc.klist_frac,klist1,klist2] = H_hr_arc.kmesh3D(kmesh,kfermiarc,'fermiarc');
switch fin_dir
case 1
klist1 =klist1(:,2);
klist2 =klist2(:,3);
case 2
klist1 =klist1(:,1);
klist2 =klist2(:,3);
case 3
klist1 =klist1(:,1);
klist2 =klist2(:,2);
end
[H00_H01_cell_list_1,H00_H01_cell_list_2] = H_hr_arc.H00_H11_cell_list_gen(fin_dir,principle_layer);
GREENCAR = HR.GREENCAR_gen(w_arc,eta,H00_H01_cell_list_1,H00_H01_cell_list_2,mode);
DOSCAR_b = HR.DOSCAR_gen(GREENCAR.bulk,'green');
DOSCAR_l = HR.DOSCAR_gen(GREENCAR.surf_l,'green');
DOSCAR_r = HR.DOSCAR_gen(GREENCAR.surf_r,'green');
DOSCAR_b = reshape(DOSCAR_b,kmesh(2),kmesh(1));
DOSCAR_l = reshape(DOSCAR_l,kmesh(2),kmesh(1));
DOSCAR_r = reshape(DOSCAR_r,kmesh(2),kmesh(1));
end
