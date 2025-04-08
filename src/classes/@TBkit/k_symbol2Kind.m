function Kind = k_symbol2Kind(H_hk,k_symbol)
k_symbol = string(k_symbol);
symvar_list = symvar(str2sym(k_symbol));
for i = 1:length(symvar_list)
str_tmp = string(symvar_list(i));
switch str_tmp
case {'k_x','k_X','K_X','K_x','kx','kX','KX','Kx','x','X'}
k_symbol = strrep(k_symbol,str_tmp,'k_x');
case {'k_y','k_Y','K_y','K_Y','ky','kY','Ky','KY','y','Y'}
k_symbol = strrep(k_symbol,str_tmp,'k_y');
case {'k_z','k_Z','K_z','K_Z','kz','kZ','Kz','KZ','z','Z'}
k_symbol = strrep(k_symbol,str_tmp,'k_z');
case {'k_w','k_W','K_w','K_W','kw','kW','Kw','KW'}
k_symbol_str = strrep(k_symbol_str,str_tmp,'k_w');
end
end
coeff_tmp = coeffs(str2sym(k_symbol));
str_2_compare = string((str2sym(k_symbol)/coeff_tmp));
Kind = find(string(H_hk.HsymL) == str_2_compare);
end
