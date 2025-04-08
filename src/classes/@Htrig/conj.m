function H_htrig = conj(H_htrig)
[num_label,coe_label] = H_htrig.NumOrCoe();
for i =1:H_htrig.Kinds
if num_label
H_htrig.HnumL(:,:,i) = conj(H_htrig.HnumL(:,:,i));
end
if coe_label
H_htrig.HcoeL(:,:,i) = conj(H_htrig.HcoeL(:,:,i)');
end
end
H_htrig.HsymL_trig = conj(H_htrig.HsymL_trig);
end
