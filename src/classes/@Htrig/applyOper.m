function H_htrig = applyOper(H_htrig,SymOper,options)
arguments
H_htrig Htrig;
SymOper Oper = Oper();
options.fast = true;
end
[~,coe_label] = H_htrig.NumOrCoe();
if ~coe_label
H_htrig = H_htrig.init();
H_htrig = H_htrig.hermitize();
end
if length(SymOper) == 1
if ~SymOper.conjugate && ~SymOper.antisymmetry && isequal(SymOper.R,eye(3))
return;
end
[H_htrig_bk,H_htrig]  = H_htrig.applyR(inv(SymOper.R));
if isnan(SymOper.U)
end
H_htrig_bk  = H_htrig_bk.applyU(SymOper.U,SymOper.conjugate,SymOper.antisymmetry);
Equationlist_r = (real(H_htrig.HcoeL - H_htrig_bk.HcoeL) == 0);
Equationlist_i = (imag(H_htrig.HcoeL - H_htrig_bk.HcoeL) == 0);
Equationlist_r = Htrig.isolateAll(Equationlist_r);
Equationlist_i = Htrig.isolateAll(Equationlist_i);
HcoeLtmp = H_htrig.HcoeL ;
HcoeLtmp_r = subs(real(HcoeLtmp),lhs(Equationlist_r),rhs(Equationlist_r));
HcoeLtmp_i = subs(imag(HcoeLtmp),lhs(Equationlist_i),rhs(Equationlist_i));
H_htrig.HcoeL = HcoeLtmp_r + 1i*HcoeLtmp_i;
else
for i = 1:length(SymOper)
H_htrig = H_htrig.applyOper(SymOper(i));
end
end
end
