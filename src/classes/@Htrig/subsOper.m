function H_htrig_bk = subsOper(H_htrig,SymOper)
arguments
H_htrig Htrig;
SymOper Oper = Oper();
end
if isequal(sym(zeros(size(H_htrig.HcoeL))),H_htrig.HcoeL)
coe_label = false;
else
coe_label = true;
end
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
end
end
