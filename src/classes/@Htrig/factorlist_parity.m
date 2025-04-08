function Factorlist_parity = factorlist_parity(H_htrig)
syms k_x k_y k_z real;
HsymC = H_htrig.HsymL_trig;
HsymC = subs(HsymC,[k_x k_y k_z],-[k_x k_y k_z]);
Factorlist_parity = simplify(HsymC./H_htrig.HsymL_trig);
end
