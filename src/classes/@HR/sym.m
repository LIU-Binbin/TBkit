function Hsym = sym(H_hr,options)
% SYM Generate symbolic Hamiltonian representation
%
%   HSYM = SYM(H_HR,OPTIONS) creates symbolic Hamiltonian
%
%   Inputs:
%       H_hr - HR object
%       options.cartesian - Use cartesian coordinates [default: true]
%       options.simple - Skip tij multiplication [default: false]
%   Output:
%       Hsym - Symbolic Hamiltonian matrix
%
%   Notes:
%       - Converts to full matrix if sparse
%       - Handles both cartesian and fractional coordinates
%       - Includes orbital symmetry factors
arguments
H_hr HR;
options.cartesian = true;
options.simple = false;
end
if strcmp(H_hr.Type,'list')
H_hr = H_hr.rewind();
elseif strcmp(H_hr.Type,'sparse')
H_hr = H_hr.full();
end
Hsym = zeros(H_hr.WAN_NUM,'sym');
H_hr = H_hr.tjmti_gen('sym');
vectorList = double(H_hr.vectorL);
[H_hr.num,~] = H_hr.NumOrCoe;
if H_hr.num
H_hr.HcoeL = sym(H_hr.HnumL);
end
if options.cartesian
tij_mat_k = H_hr.tjmti{3};
VarsUsing = H_hr.VarsSeqLcart(1:H_hr.Dim);
exp_pre = exp(1i...
*VarsUsing*(...
vectorList*H_hr.Rm)');
for i = 1:H_hr.NRPTS
Hsym = Hsym+H_hr.HcoeL(:,:,i)*exp_pre(i);
end
if options.simple
else
Hsym = Hsym.*tij_mat_k;
end
else
tij_mat_s = H_hr.tjmti{4};
VarsUsing = H_hr.VarsSeqLfrac(1:H_hr.Dim);
exp_pre = exp(1i*2*pi...
*VarsUsing*(...
vectorList)');
for i = 1:H_hr.NRPTS
Hsym = Hsym+H_hr.HcoeL(:,:,i)*exp_pre(i);
end
if options.simple
else
Hsym = Hsym.*tij_mat_s;
end
end
end
