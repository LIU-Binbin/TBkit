function [H_hk_R,H_hk] = applyRU(H_hk,SymOper)
if ~SymOper.conjugate && ~SymOper.antisymmetry && isequal(SymOper.R,eye(3))
H_hk_R = H_hk;
return;
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
if SymOper.conjugate
R_inv = -inv(SymOper.R);
else
R_inv = inv(SymOper.R);
end
H_hk_R = H_hk;
matcell = H_hk.matgen(R_inv);
if num_label
for i = 0:H_hk.Degree
Orderlist = HK.orderlist(i);
H_hk_R.HnumL(:,:,Orderlist) = HK.matrixtimespage(matcell{i+1}.',H_hk_R.HnumL(:,:,Orderlist));
end
end
if coe_label
for i = 0:H_hk.Degree
Orderlist = HK.orderlist(i);
H_hk_R.HcoeL(:,:,Orderlist) = HK.matrixtimespage(matcell{i+1}.',H_hk_R.HcoeL(:,:,Orderlist));
end
end
if isnan(SymOper.U)
else
U = SymOper.U;
end
U_inv = inv(U);
if coe_label == true
HcoeLtmp = H_hk_R.HcoeL;
if SymOper.conjugate
HcoeLtmp = conj(HcoeLtmp);
end
if SymOper.antisymmetry
HcoeLtmp = -HcoeLtmp;
end
HcoeLtmp = HK.page_mtimes_matrix(HK.matrix_mtimes_page(U,HcoeLtmp),U_inv);
H_hk_R.HcoeL = (HcoeLtmp);
end
if num_label == true
U_page = repmat(U,[1 1 size(H_hk_R.HnumL,3)]);
U_inv_page = repmat(U_inv,[1 1 size(H_hk_R.HnumL,3)]);
HnumLtmp = H_hk_R.HnumL;
if conjugate
HnumLtmp = conj(HnumLtmp);
HnumLtmp = HK.matrixtimepage(HK.factorlist_parity(H_hk_R.Degree),HnumLtmp);
end
if antisymmetry
HnumLtmp = -HnumLtmp;
end
HnumLtmp = pagemtimes(pagemtimes(U_page,HnumLtmp),U_inv_page);
H_hk_R.HnumL = HnumLtmp;
end
end
