function H_hr = set_overlap(H_hr,amp,si,sj,vector_list,mode)
% SET_OVERLAP Set overlap terms in HR object
%
%   H_HR = SET_OVERLAP(H_HR,AMP,SI,SJ,VECTOR_LIST,MODE) sets overlap terms
%
%   Inputs:
%       H_hr - HR object to modify
%       amp - Overlap amplitude(s)
%       si - Initial orbital index(es)
%       sj - Final orbital index(es)
%       vector_list - List of lattice vectors
%       mode - Operation mode ('set', 'add', 'sym', 'symadd') [default: 'set']
%   Output:
%       H_hr - Modified HR object
%
%   Notes:
%       - High-level interface for setting overlaps
%       - Handles scalar and matrix inputs
%       - Dispatches to set_overlap_mat or set_overlap_single as needed
%       - Supports WAN_NUM/2 subblock operations
if nargin <6
mode = 'set';
end
WAN_NUM_half = H_hr.WAN_NUM/2;
[n_vector,~] = size(vector_list);
for i = 1:n_vector
vector = vector_list(i,:);
if length(si) == length(sj) && length(si)>1 && strcmp(H_hr.Type,'mat')
WANNUM = H_hr.WAN_NUM;
switch mode
case {'set','add'}
amp_mat =zeros(WANNUM);
case {'sym','symadd'}
amp_mat =sym(zeros(WANNUM));
end
amp_mat(sub2ind([WANNUM,WANNUM],si,sj)) = amp;
H_hr = H_hr.set_overlap_mat(amp_mat,vector,mode);
elseif length(si) == length(sj) && length(si)>1 && strcmp(H_hr.Type,'list')
for j = 1:length(si)
H_hr = H_hr.set_overlap_single(amp(j),si(j),sj(j),vector,mode);
end
else
length_amp = length(amp);
switch length_amp
case 1
H_hr = H_hr.set_overlap_single(amp,si,sj,vector,mode);
case 2
H_hr = H_hr.set_overlap_single(amp(1,1),si,sj,vector,mode);
H_hr = H_hr.set_overlap_single(amp(1,2),si,sj+WAN_NUM_half,vector,mode);
H_hr = H_hr.set_overlap_single(amp(2,1),si+WAN_NUM_half,sj,vector,mode);
H_hr = H_hr.set_overlap_single(amp(2,2),si+WAN_NUM_half,sj+WAN_NUM_half,vector,mode);
case 4
disp('if needed');
case 8
disp('if needed');
otherwise
H_hr = H_hr.set_overlap_mat(amp,vector,mode);
end
end
end
end
