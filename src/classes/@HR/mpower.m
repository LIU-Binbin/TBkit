function H_hr = mpower(A,B)
if isa(A,'HR') && isa(B,'numeric')
if length(B) >1
error('only support a number');
end
H_hr = A;
c = B;
for i = 1:H_hr.NRPTS
H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i)^c;
H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i)^c;
end
elseif isa(B,'HR') && isa(A,'numeric')
if length(A) >1
error('only support a number');
end
H_hr = B;
c = A;
for i = 1:H_hr.NRPTS
H_hr.HcoeL(:,:,i) = c^H_hr.HcoeL(:,:,i);
H_hr.HnumL(:,:,i) = c^H_hr.HnumL(:,:,i);
end
else
error('wrong input');
end
end
