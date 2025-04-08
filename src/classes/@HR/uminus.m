function H_hr = uminus(H_hr)
for i = 1:H_hr.NRPTS
H_hr.HcoeL(:,:,i) = -H_hr.HcoeL(:,:,i);
H_hr.HnumL(:,:,i) = -H_hr.HnumL(:,:,i);
end
end
