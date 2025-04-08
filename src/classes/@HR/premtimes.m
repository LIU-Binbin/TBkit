function H_hr = premtimes(A,B)
if isa(A,'HR') && ~isa(B,'HR')
H_hr1 = A;
H_hr = H_hr1;
if isa(B,'sym')
for i = 1:H_hr.NRPTS
H_hr.HcoeL(:,:,i) = B*H_hr.HcoeL(:,:,i);
end
elseif isa(B,'numeric')
for i = 1:H_hr.NRPTS
H_hr.HnumL(:,:,i) = B*H_hr.HnumL(:,:,i) ;
end
else
end
end
end
