function sym_coe_list = coeff_extract(sym_term)
str_list_tmp = strsplit(string(sym_term),'*');
sym_coe = sym(1);
for i = 1:length(str_list_tmp)
if Htrig.strcontain(str_list_tmp{i},['x','y','z'])
sym_coe_list(i) = sym_coe * str2sym(str_list_tmp{i});
end
end
end
