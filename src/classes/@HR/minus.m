function H_hr = minus(A,B)
if isa(A,'HR') && isa(B,'HR')
H_hr1 = A;
H_hr2 = B;
if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
error('WAN_NUM different');
end
H_hr = H_hr1;
for i = 1:H_hr2.NRPTS
vector = H_hr2.vectorL(i,:);
amp = -H_hr2.HnumL(:,:,i);
amp_sym = -H_hr2.HcoeL(:,:,i);
H_hr = H_hr.set_hop_mat(amp,vector,'add');
H_hr = H_hr.set_hop_mat(amp_sym,vector,'symadd');
end
elseif isa(A,'HR') && ~isa(B,'HR')
H_hr1 = A;
H_hr = H_hr1;
i = H_hr.Line_000;
if isa(B,'sym')
H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i) - B;
elseif isa(B,'numeric')
H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i) - B;
else
end
elseif ~isa(A,'HR') && isa(B,'HR')
H_hr1 = B;
H_hr = H_hr1;
if isa(A,'sym')
for i = 1:H_hr.NRPTS
H_hr.HcoeL(:,:,i) = A + -H_hr.HcoeL(:,:,i) ;
end
elseif isa(A,'numeric')
for i = 1:H_hr.NRPTS
H_hr.HnumL(:,:,i) = A + -H_hr.HnumL(:,:,i) ;
end
else
end
end
end
