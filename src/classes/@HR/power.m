function H_hr = power(H_hr,b)
for i = 1:H_hr.NRPTS
H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i).^b;
H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i).^b;
end
end
