function H_hr = ctranspose(H_hr)
H_hr = H_hr.dualize();
if H_hr.vectorhopping
H_hr.AvectorL = H_hr.AvectorL(H_hr.Duality_vector_dist,:);
H_hr.BvectorL = -H_hr.BvectorL(H_hr.Duality_vector_dist,:);
CLr = H_hr.CvectorL(1:end/2,:);
CLi = H_hr.CvectorL(end/2+1:end,:);
H_hr.CvectorL = [CLr(H_hr.Duality_vector_dist,:);-CLi(H_hr.Duality_vector_dist,:)];
else
if strcmp(H_hr.Type,'list')
if H_hr.num
H_hr.HnumL = conj(H_hr.HnumL(H_hr.Duality_vector_dist));
end
if H_hr.coe
H_hr.HcoeL = conj(H_hr.HcoeL(H_hr.Duality_vector_dist));
end
else
if H_hr.coe
H_hr.HcoeL = pagectranspose(H_hr.HcoeL(:,:,H_hr.Duality_vector_dist));
end
if H_hr.num
H_hr.HnumL = pagectranspose(H_hr.HnumL(:,:,H_hr.Duality_vector_dist));
end
end
end
end
