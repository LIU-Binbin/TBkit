function Factorlist_parity = factorlist_parity(H_hr)
syms k_x k_y k_z real;
HsymC = H_hr.HsymL_trig;
HsymC = subs(HsymC,[k_x k_y k_z],-[k_x k_y k_z]);
Factorlist_parity = simplify(HsymC./H_hr.HsymL_trig);
end
