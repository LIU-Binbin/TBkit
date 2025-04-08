function [sym_coe,sym_single,str_single] = coeff_extract(sym_term,Dim)
arguments
sym_term
Dim = 3;
end
VarsSeqLcart = [sym('k_x'),sym('k_y'),sym('k_z'),sym('k_w')];
[sym_coe,sym_single] = coeffs(sym_term,VarsSeqLcart(1:Dim));
str_single = string(sym_single);
end
