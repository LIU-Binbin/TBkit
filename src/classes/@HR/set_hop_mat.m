function H_hr = set_hop_mat(H_hr,amp,vector,mode)
% SET_HOP_MAT Set hopping matrix elements in HR object
%
%   H_HR = SET_HOP_MAT(H_HR,AMP,VECTOR,MODE) sets hopping matrix elements
%
%   Inputs:
%       H_hr - HR object to modify
%       amp - Hopping amplitude matrix
%       vector - Lattice vector for hopping
%       mode - Operation mode ('set', 'add', 'sym', 'symadd')
%   Output:
%       H_hr - Modified HR object
%
%   Notes:
%       - Handles all storage formats ('sparse', 'list', 'mat')
%       - Supports both numeric and symbolic operations
%       - Automatically adds new vectors if needed
V = H_hr.vectorL;
switch H_hr.Type
case 'list'
sizeamp = size(amp);
for n = 1:numel(amp)
[hi,hj] = ind2sub(sizeamp,n);
if ~isequal(sym(amp(n)),sym(0))
H_hr = set_hop_single(H_hr,amp(n),hi,hj,vector,mode);
end
end
case 'sparse'
[~,seq]=ismember((vector),V,'rows');
if seq == 0
seq = H_hr.NRPTS +1;
H_hr = H_hr.add_empty_one(vector);
end
switch mode
case 'set'
%    fprintf('%d %d %d %d %d \n',vector(1),vector(2),vector(3));
H_hr.HnumL{seq} = amp ;
case 'add'
H_hr.HnumL{seq} = H_hr.HnumL{seq} + amp  ;
case 'sym'
error('not be implemented yet');
case 'symadd'
error('not be implemented yet');
end
otherwise
[~,seq]=ismember((vector),V,'rows');
if seq == 0
seq = H_hr.NRPTS +1;
H_hr = H_hr.add_empty_one(vector);
end
switch mode
case 'set'
zeros_num_mat = zeros(H_hr.WAN_NUM);
if H_hr.HnumL(:,:,seq) ~= zeros_num_mat
warning('May be you should use add mode on this NRPT');
Format = fold(@strcat,repmat("
Input = num2cell(vector(1,:));
fprintf(strcat(Format," \n"),Input{:});
end
H_hr.HnumL(:,:,seq) = amp ;
case 'add'
H_hr.HnumL(:,:,seq) = H_hr.HnumL(:,:,seq) + amp  ;
case 'sym'
zeros_coe_mat = sym(zeros(H_hr.WAN_NUM));
if H_hr.HcoeL(:,:,seq) ~= zeros_coe_mat
warning('May be you should use symadd mode on this NRPT ');
Format = fold(@strcat,repmat("
Input = num2cell(vector(1,:));
fprintf(strcat(Format," \n"),Input{:});
end
H_hr.HcoeL(:,:,seq) =  amp  ;
case 'symadd'
H_hr.HcoeL(:,:,seq) = H_hr.HcoeL(:,:,seq) + amp  ;
end
end
end
