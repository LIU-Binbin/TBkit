function H_hr = H_TB_gen_SK(H_hr,options)
% H_TB_GEN_SK Generate tight-binding Hamiltonian (Slater-Koster)
%
%   H_hr = H_TB_GEN_SK(H_hr,options) constructs a tight-binding
%   Hamiltonian using Slater-Koster parameters.
%
%   INPUT ARGUMENTS:
%       H_hr - HR object with neighbor information
%       options - Structure with parameters:
%           level_cut: Maximum hopping level
%           onsite: Include onsite terms (logical)
%           per_dir: Periodic directions
%           chiral: Chiral filtering (logical)
%           spin: Spin filtering (logical)
%
%   OUTPUT ARGUMENTS:
%       H_hr - HR object with generated Hamiltonian
%
%   NOTES:
%       - Requires precomputed nn_store_smart
%       - Uses symbolic representation for coefficients
%       - Includes progress reporting
%
%   SEE ALSO:
%       HR, set_hop
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

arguments
H_hr HR;
options.level_cut {mustBeInteger} = 1;
options.onsite logical = false;
options.per_dir double = [1,1,1];
options.chiral logical = false;
options.spin logical = false;
options.method  = 'nn_sk_smart';
options.rough = false;
options.vectorL = [];
end
if isempty(H_hr.nn_store_smart)
if isempty(H_hr.nn_store_smart)
error('You have not run the function: nn_sk_smart');
else
if strcmp(options.method,'nn_sparse')
error('You should use the key value: ''method'',''nn_smart'' as input.');
end
end
end
N_orbit = H_hr.WAN_NUM;
N_tot = H_hr.WAN_NUM;
H_hr.HcoeL(:,:,1)  = sym(zeros(N_orbit));
H_hr.HnumL(:,:,1)  = zeros(N_orbit);
level_cut = options.level_cut;
per_dir = options.per_dir;
if options.onsite
for i=1:N_orbit
onsite_sym_name = "E__"+string(H_hr.elementL(i,1))+"_"...
+string(H_hr.quantumL(i,2));...
H_hr = H_hr.set_hop_single(sym(onsite_sym_name,'real'),i,i,[0 0 0],'sym');
end
end
pb = CmdLineProgressBar('Setting ');
for j=1:N_orbit
pb.print(j,N_orbit,' th orbital ...');
spin1 = H_hr.quantumL(j,4);
element1 =  H_hr.elementL(j);
for i=1:N_tot
spin2 = H_hr.quantumL(i,4);
element2 =  H_hr.elementL(i);
set_or_not = true;
if options.chiral
if element1 == element2
set_or_not = false;
end
end
if options.spin
if spin1 == spin2
set_or_not = false;
end
end
if set_or_not
nn = H_hr.nn_store_smart(i,j).nn;
for k = 1:size(nn,1)
if nn(k).nn_level <= level_cut
vector_init = nn(k).R_vector.*per_dir;
H_hr = H_hr.set_hop(nn(k).hop,i,j,vector_init,'sym');
if H_hr.overlap
H_hr = H_hr.set_hop(nn(k).overlap,i,j,vector_init,'sym');
end
end
end
end
end
end
pb.delete();
H_hr.coe = true;
H_hr.num = false;
end
