function [klist_cart,klist_frac,gap_list,fig] = findnodes2(H_hr,kzf_list,Noccupy,tolerance)
% FINDNODES2 Locate band crossing nodes along kz direction
%
%   [klist_cart,klist_frac,gap_list,fig] = FINDNODES2(H_hr,kzf_list,Noccupy,tolerance)
%   searches for band crossings by scanning kz values.
%
%   INPUT ARGUMENTS:
%       H_hr - Hamiltonian in HR format
%       kzf_list - List of kz values to scan (default: 0:0.01:1)
%       Noccupy - Number of occupied bands (default: WAN_NUM/2)
%       tolerance - Gap threshold for node detection (default: 0.01)
%
%   OUTPUT ARGUMENTS:
%       klist_cart - Node positions in cartesian coordinates
%       klist_frac - Node positions in fractional coordinates
%       gap_list - Band gaps at node positions
%       fig - Figure handle for band gap plot
%
%   NOTES:
%       - Plots band gap vs kz during search
%       - Returns only points where gap < tolerance
%
%   SEE ALSO:
%       HR, EIGENCAR_gen
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]
if nargin < 4
tolerance = 0.01;
end
if nargin <3
Noccupy = H_hr.WAN_NUM/2;
end
if nargin <2
kzf_list = 0:0.01:1;
end
import TBkit_tool.*;
[fig,ax] = creat_figure();
xlabel(ax,'kz_f');
ylabel(ax,'min gap in k_z plane');
title(ax,'findnode k_z plane mode');
Gk_ = (eye(3)*2*pi)/H_hr.Rm;
klist_cart = [];
klist_frac = [];
gap_list =[];
for i = 1:length(kzf_list)
H_hr.klist_frac(:,3) = kzf_list(i);
EIGENCAR = H_hr.EIGENCAR_gen();
[bandgap,label] = min(EIGENCAR(Noccupy+1,:)- EIGENCAR(Noccupy,:));
scatter(ax,kzf_list(i),EIGENCAR(Noccupy,label));
scatter(ax,kzf_list(i),EIGENCAR(Noccupy+1,label));
drawnow limitrate;
if bandgap < tolerance
fprintf('find it: %7.4f eV (%6.3f, %6.3f,%6.3f)\n',bandgap,H_hr.klist_frac(label,1),H_hr.klist_frac(label,2),kzf_list(i));
klist_cart = [ klist_cart ;H_hr.klist_cart(label,1),H_hr.klist_cart(label,2),kzf_list(i)*(Gk_(3,3))];
klist_frac = [klist_frac ;H_hr.klist_frac(label,1),H_hr.klist_frac(label,2),kzf_list(i)];
gap_list =[gap_list;bandgap];
else
fprintf('min gap: %7.4f eV (%6.3f, %6.3f,%6.3f)\n',bandgap,H_hr.klist_frac(label,1),H_hr.klist_frac(label,2),kzf_list(i));
end
end
end
