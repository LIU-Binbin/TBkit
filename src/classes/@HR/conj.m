function H_hr = conj(H_hr)
H_hr = H_hr.dualize();
H_hr.HcoeL = conj(H_hr.HcoeL(:,:,H_hr.Duality_vector_dist));
H_hr.HnumL = conj(H_hr.HnumL(:,:,H_hr.Duality_vector_dist));
end
