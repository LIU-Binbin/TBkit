function H_hk = Degree2Kinds(H_hk)
VarUsing = H_hk.VarsSeqLcart(1:H_hk.Dim);
str_tmp = string(expand((1+fold(@plus,VarUsing))^H_hk.Degree));
str_cell = strsplit(str_tmp,'+');
H_hk.Kinds = length(str_cell);
coeff_list = zeros(H_hk.Kinds,1);
Degree_list = zeros(H_hk.Kinds,1);
for i = 1:H_hk.Kinds
tmpsym = str2sym(str_cell{i});
Degree_list(i) = polynomialDegree(tmpsym);
coeff_list(i) = coeffs(tmpsym);
end
[~,sort_label] = sort(Degree_list);
coeff_list = coeff_list(sort_label);
str_cell = str_cell(sort_label);
for i = 1:H_hk.Kinds
H_hk.HsymL(i) = str2sym(str_cell{i})/coeff_list(i);
end
H_hk.HstrL = string(H_hk.HsymL);
temp_strL = H_hk.HstrL;
temp_strL = strrep(temp_strL,'k_x','x');
temp_strL = strrep(temp_strL,'k_y','y');
temp_strL = strrep(temp_strL,'k_z','z');
for i = 1:H_hk.Kinds
H_hk.HstrL(i) = HK.quadk2K(temp_strL(i));
end
H_hk.HsymL_xyz = str2sym(temp_strL);
end
