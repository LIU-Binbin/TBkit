function H_hk = applyU(H_hk,U,conjugate ,antisymmetry )
arguments
H_hk HK;
U = nan;
conjugate logical =false;
antisymmetry logical = false;
end
if isa(U,'Oper')
conjugate = U.conjugate;
antisymmetry = U.antisymmetry;
U = U.U;
end
if isequal(zeros(size(H_hk.HnumL)),H_hk.HnumL)
num_label = false;
else
num_label = true;
end
if isequal(sym(zeros(size(H_hk.HcoeL))),H_hk.HcoeL)
coe_label = false;
else
coe_label = true;
end
U_inv = inv(U);
if coe_label == true
HcoeLtmp = H_hk.HcoeL;
if conjugate
HcoeLtmp = conj(HcoeLtmp);
HcoeLtmp = HK.matrixtimespage(HK.factorlist_parity(H_hk.Degree),HcoeLtmp);
end
if antisymmetry
HcoeLtmp = -HcoeLtmp;
end
HcoeLtmp = HK.page_mtimes_matrix(HK.matrix_mtimes_page(U,HcoeLtmp),U_inv);
H_hk.HcoeL = HcoeLtmp;
end
if num_label == true
U_page = repmat(U,[1 1 size(H_hk.HnumL,3)]);
U_inv_page = repmat(U_inv,[1 1 size(H_hk.HnumL,3)]);
HnumLtmp = H_hk.HnumL;
if conjugate
HnumLtmp = conj(HnumLtmp);
HnumLtmp = HK.matrixtimepage(HK.factorlist_parity(H_hk.Degree),HnumLtmp);
end
if antisymmetry
HnumLtmp = -HnumLtmp;
end
HnumLtmp = pagemtimes(pagemtimes(U_page,HnumLtmp),U_inv_page);
H_hk.HnumL = HnumLtmp;
end
end
