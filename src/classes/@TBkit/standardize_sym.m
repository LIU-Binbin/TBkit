function sym_term = standardize_sym(sym_term)
str_tmp = string(sym_term);
str_tmp = replace(str_tmp,{'k_x','k_X','K_X','K_x','kx','kX','KX','Kx'},'k_x');
str_tmp = replace(str_tmp,{'k_y','k_Y','K_y','K_Y','ky','kY','Ky','KY'},'k_y');
str_tmp = replace(str_tmp,{'k_z','k_Z','K_z','K_Z','kz','kZ','Kz','KZ'},'k_z');
sym_term = str2sym(str_tmp);
end
