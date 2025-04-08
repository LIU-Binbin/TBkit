function H_htrig = applyU(H_htrig,U,conjugate ,antisymmetry )
arguments
H_htrig Htrig;
U =nan;
conjugate logical =false;
antisymmetry logical = false;
end
if isa(U,'Oper')
conjugate = U.conjugate;
antisymmetry = U.antisymmetry;
U = U.U;
end
[num_label,coe_label] = H_htrig.NumOrCoe();
U_inv = inv(U);
if coe_label == true
HcoeLtmp = H_htrig.HcoeL;
if conjugate
if strcmp(H_htrig.Type,'sincos')
H_htrig = H_htrig.applyR(diag([-1,-1,-1]));
HcoeLtmp = H_htrig.HcoeL;
end
HcoeLtmp = conj(HcoeLtmp);
end
if antisymmetry
HcoeLtmp = -HcoeLtmp;
end
HcoeLtmp = Htrig.page_mtimes_matrix(Htrig.matrix_mtimes_page(U,HcoeLtmp),U_inv);
H_htrig.HcoeL = HcoeLtmp;
end
if num_label == true
U_page = repmat(U,[1 1 size(H_htrig.HnumL,3)]);
U_inv_page = repmat(U_inv,[1 1 size(H_htrig.HnumL,3)]);
HnumLtmp = H_htrig.HnumL;
if conjugate
HnumLtmp = conj(HnumLtmp);
HnumLtmp = Htrig.matrixtimespage(H_htrig.factorlist_parity(),HnumLtmp);
end
if antisymmetry
HnumLtmp = -HnumLtmp;
end
HnumLtmp = pagemtimes(pagemtimes(U_page,HnumLtmp),U_inv_page);
H_htrig.HnumL = HnumLtmp;
end
end
