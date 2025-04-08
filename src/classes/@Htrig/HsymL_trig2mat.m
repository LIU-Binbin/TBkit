function mat = HsymL_trig2mat(H_htrig,HsymL_trig)
NSLAB = (H_htrig.Nslab ==0)+H_htrig.Nslab;
NS = prod(NSLAB);
Delta_ijk = Htrig.Delta_Oper(Htrig.coeff_extract(HsymL_trig));
tmpmat = meshgrid((1:NS).',(1:NS).');
NS1 = reshape(tmpmat,NS^2,1);
NS2 = reshape(tmpmat.',NS^2,1);
[i1,i2,i3] = ind2sub(NSLAB,NS1);
[j1,j2,j3] = ind2sub(NSLAB,NS2);
mat = double(reshape(Delta_ijk(i1,i2,i3,j1,j2,j3),NS,NS)).';
end
