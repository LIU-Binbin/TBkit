function H_hr = H_TB_gen_SK_sparse(H_hr,level_cut,para_filename,onsite_mode)

% H_TB_GEN_SK_SPARSE Generate sparse tight-binding Hamiltonian (Slater-Koster)
%
%   H_hr = H_TB_GEN_SK_SPARSE(H_hr,options) constructs a sparse
%   tight-binding Hamiltonian using Slater-Koster parameters.
%
%   INPUT ARGUMENTS:
%       H_hr - HR object with neighbor information
%       options - Structure with parameters:
%           level_cut: Maximum hopping level
%           onsite: Include onsite terms (logical)
%           para: Parameter structure
%           deltarule: Delta rule type (0-2)
%           alpharule: Alpha rule type (0-2)
%
%   OUTPUT ARGUMENTS:
%       H_hr - HR object with generated Hamiltonian
%
%   NOTES:
%       - Requires precomputed nn_store_smart
%       - Supports chiral and spin filtering
%       - Uses TBSK_Var_gen_sparse for parameter generation
%
%   SEE ALSO:
%       HR, TBSK_Var_gen_sparse
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

if nargin <4
onsite_mode = 0;
end
if nargin <3
para_filename = 'para.mat';
end
if nargin <2
level_cut = 1;
end
if exist(para_filename,'file')
load(para_filename);
else
error('in sparse mode, we give value for the parameter first!');
end
N_orbit = H_hr.WAN_NUM;
nn_sparse_num = zeros(size(H_hr.nn_store_smart,1),5);
H_hr.vectorL = unique(H_hr.nn_store_smart(:,3:5),'rows');
Nhopping = size(H_hr.nn_store_smart,1);
for i = 1: Nhopping
nn_sparse_num(i,1) = H_hr.nn_store_smart(i,1);
nn_sparse_num(i,2) = H_hr.nn_store_smart(i,2);
[~,nn_sparse_num(i,3)] = ismember(...
[H_hr.nn_store_smart(i,3),...
H_hr.nn_store_smart(i,4),...
H_hr.nn_store_smart(i,5)],...
H_hr.vectorL,'rows');
orb1  = H_hr.sites(nn_sparse_num(i,1)).orb;
orb2  = H_hr.sites(nn_sparse_num(i,2)).orb;
if strcmp(orb1,'p')
if strcmp(orb2,'s')
orb1 = 's';
orb2 = 'p';
end
end
nn_sparse_num(i,4) = H_hr.nn_store_smart(i,7);
if nn_sparse_num(i,4)  <= level_cut
nn_sparse_num(i,5)=nn_sparse_num(i,5)+...
H_hr.nn_store_smart(i,8)*subs(str2sym("V"+orb1+orb2+'S_'+string(H_hr.nn_store_smart(i,7))))+...
H_hr.nn_store_smart(i,9)*subs(str2sym("V"+orb1+orb2+'P_'+string(H_hr.nn_store_smart(i,7))));...
else
nn_sparse_num(i,5)=0;
end
end
H_hr.HnumL = [];
H_hr.HnumL{H_hr.NRPTS}  = sparse(N_orbit,N_orbit);
for i =1:H_hr.NRPTS -1
H_hr.HnumL{i}  = sparse(N_orbit,N_orbit);
end
if onsite_mode == 1
for j=1:N_orbit
onsite_sym_name = "E_onsite_"+string(H_hr.quantumL(j,1))+"_"...
+string(H_hr.quantumL(j,2));...
H_hr.HnumL{i}(j,j) = subs(str2sym(onsite_sym_name));
end
end
for Nh=1: Nhopping
fprintf("setting (%4d/%4d) th hopping ... \n",Nh, Nhopping);
if nn_sparse_num(Nh,4)  <= level_cut
H_hr.HnumL{nn_sparse_num(Nh,3)}(nn_sparse_num(Nh,1),nn_sparse_num(Nh,2)) = ...
nn_sparse_num(Nh,5);
end
end
H_hr.nn_sparse_n = nn_sparse_num;
end
