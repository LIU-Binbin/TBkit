function H_htrig = hermitize(H_htrig)
[num_label,coe_label] = H_htrig.NumOrCoe();
H_htrig_bk = H_htrig';
if coe_label
Equationlist_r = real(H_htrig.HcoeL - H_htrig_bk.HcoeL) == 0;
Equationlist_i = imag(H_htrig.HcoeL - H_htrig_bk.HcoeL) == 0;
Equationlist_r = Htrig.isolateAll(Equationlist_r);
Equationlist_i = Htrig.isolateAll(Equationlist_i);
HcoeLtmp = subs(H_htrig.HcoeL,lhs(Equationlist_r),rhs(Equationlist_r));
HcoeLtmp = subs(HcoeLtmp,lhs(Equationlist_i),rhs(Equationlist_i));
H_htrig.HcoeL = HcoeLtmp;
end
if num_label
H_htrig.HnumL = (H_htrig_bk.HnumL + H_htrig.HnumL )/2;
end
end
