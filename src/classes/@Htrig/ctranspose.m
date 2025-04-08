function H_htrig = ctranspose(H_htrig)
[num_label,coe_label] = H_htrig.NumOrCoe();
for i =1:H_htrig.Kinds
if num_label
H_htrig.HnumL(:,:,i) = H_htrig.HnumL(:,:,i)';
end
if coe_label
H_htrig.HcoeL(:,:,i) = H_htrig.HcoeL(:,:,i)';
end
end
if strcmp(H_htrig.Type,'exp')
H_htrig =H_htrig.dualize();
if num_label
H_htrig.HnumL(:,:,:) = H_htrig.HnumL(:,:,H_htrig.Duality_vector_dist);
end
if coe_label
H_htrig.HcoeL(:,:,:) = H_htrig.HcoeL(:,:,H_htrig.Duality_vector_dist);
end
end
end
