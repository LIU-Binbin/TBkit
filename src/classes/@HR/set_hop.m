function H_hr = set_hop(H_hr,amp,hi,hj,vector_list,mode)
% SET_HOP Set hopping terms in HR object
%
%   H_HR = SET_HOP(H_HR,AMP,HI,HJ,VECTOR_LIST,MODE) sets hopping terms
%
%   Inputs:
%       H_hr - HR object to modify
%       amp - Hopping amplitude(s)
%       hi - Initial orbital index(es)
%       hj - Final orbital index(es)
%       vector_list - List of lattice vectors
%       mode - Operation mode ('set', 'add', 'sym', 'symadd')
%   Output:
%       H_hr - Modified HR object
%
%   Notes:
%       - High-level interface for setting hoppings
%       - Handles scalar and matrix inputs
%       - Dispatches to appropriate lower-level functions
arguments
H_hr;
amp;
hi;
hj;
vector_list;
mode = 'set';
end
WAN_NUM_half = H_hr.WAN_NUM/2;
[n_vector,~] = size(vector_list);
if strcmp(mode,'sym')||strcmp(mode,'symadd')
H_hr.num = false;
H_hr.coe = true;
end
for i = 1:n_vector
vector = vector_list(i,:);
if length(hi) == length(hj) && length(hi)>1 && strcmp(H_hr.Type,'mat')
WANNUM = H_hr.WAN_NUM;
switch mode
case {'set','add'}
amp_mat =zeros(WANNUM);
case {'sym','symadd'}
amp_mat =sym(zeros(WANNUM));
end
amp_mat(sub2ind([WANNUM,WANNUM],hi,hj)) = amp;
H_hr = H_hr.set_hop_mat(amp_mat,vector,mode);
elseif  length(hi) == length(hj) && strcmp(H_hr.Type,'sparse')
WANNUM = H_hr.WAN_NUM;
amp_mat = sparse(hi,hj,amp,WANNUM,WANNUM);
H_hr = H_hr.set_hop_mat(amp_mat,vector,mode);
elseif length(hi) == length(hj) && length(hi)>1 && strcmp(H_hr.Type,'list')
for j = 1:length(hi)
H_hr = H_hr.set_hop_single(amp(j),hi(j),hj(j),vector,mode);
end
else
length_amp = length(amp);
switch length_amp
case 1
H_hr = H_hr.set_hop_single(amp,hi,hj,vector,mode);
case 2
H_hr = H_hr.set_hop_single(amp(1,1),hi,hj,vector,mode);
H_hr = H_hr.set_hop_single(amp(1,2),hi,hj+WAN_NUM_half,vector,mode);
H_hr = H_hr.set_hop_single(amp(2,1),hi+WAN_NUM_half,hj,vector,mode);
H_hr = H_hr.set_hop_single(amp(2,2),hi+WAN_NUM_half,hj+WAN_NUM_half,vector,mode);
case 4
disp('if needed');
case 8
disp('if needed');
otherwise
H_hr = H_hr.set_hop_mat(amp,vector,mode);
end
end
end
end
