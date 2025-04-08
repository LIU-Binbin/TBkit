function H_hk = ctranspose(H_hk)
for i =1:H_hk.Kinds
H_hk.HnumL(:,:,i) = H_hk.HnumL(:,:,i)';
H_hk.HcoeL(:,:,i) = H_hk.HcoeL(:,:,i)';
end
end
