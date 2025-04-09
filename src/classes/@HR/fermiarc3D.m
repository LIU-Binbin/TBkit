function [DOSCAR_l,DOSCAR_b,DOSCAR_r,klist1,klist2,E_list] = fermiarc3D(H_hr,w_range,fin_dir,kmesh,kfermiarc,options)
% FERMIARC3D Calculate 3D Fermi arc surface states
%
%   [DOSCAR_l,DOSCAR_b,DOSCAR_r,klist1,klist2,E_list] = FERMIARC3D(H_hr,w_range,fin_dir,kmesh,kfermiarc,options)
%   computes surface Green's functions and DOS over energy range.
%
%   INPUT ARGUMENTS:
%       H_hr - Bulk Hamiltonian in HR format
%       w_range - Energy range [start,end,n] or vector
%       fin_dir - Termination direction (1-3, default: 3)
%       kmesh - k-point mesh dimensions (default: [100 100])
%       kfermiarc - k-path basis (default: standard square path)
%       options - Structure with parameters:
%           principle_layer: Number of layers (default: 1)
%           eta: Broadening (default: 0.01)
%           fermi: Fermi level (default: 0)
%           mode: Calculation mode (default: 'Green_iter')
%
%   OUTPUT ARGUMENTS:
%       DOSCAR_l - Left surface DOS (3D array)
%       DOSCAR_b - Bulk DOS (3D array)
%       DOSCAR_r - Right surface DOS (3D array)
%       klist1,klist2 - k-point coordinates
%       E_list - Energy values
%
%   NOTES:
%       - Returns 3D arrays (kmesh x energy) for each DOS type
%       - Supports both linear range and explicit energy lists
%
%   SEE ALSO:
%       HR, fermiarc, GREENCAR_gen
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]
arguments
H_hr HR;
w_range double;
fin_dir{mustBeMember(fin_dir,[1,2,3])} =3;
kmesh double = [100 100];
kfermiarc = [-0.5, -0.5, 0;...
1.0 , 0,   0;...
0   , 1.0,   0;...
];
options.principle_layer = 1;
options.eta = 0.01;
options.fermi = 0;
options.mode = 'Green_iter';
end
if length(w_range) == 3
nw = w_range(3);
w_range = linspace(w_range(1),w_range(2),w_range(3))+options.fermi;
else
w_range = w_range + +options.fermi;
nw = length(w_range);
end
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
[H00_H01_cell_list_1,H00_H01_cell_list_2] = H_hr_arc.H00_H11_cell_list_gen(fin_dir,options.principle_layer);
GREENCAR = HR.GREENCAR_gen(w_range,options.eta,H00_H01_cell_list_1,H00_H01_cell_list_2,options.mode);
DOSCAR_b = HR.DOSCAR_gen(GREENCAR.bulk,'green');
DOSCAR_l = HR.DOSCAR_gen(GREENCAR.surf_l,'green');
DOSCAR_r = HR.DOSCAR_gen(GREENCAR.surf_r,'green');
DOSCAR_b = reshape(DOSCAR_b.',kmesh(2),kmesh(1),nw);
DOSCAR_l = reshape(DOSCAR_l.',kmesh(2),kmesh(1),nw);
DOSCAR_r = reshape(DOSCAR_r.',kmesh(2),kmesh(1),nw);
E_list = w_range;
end
