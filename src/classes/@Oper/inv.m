function SymOper = inv(SymOper)
if isnan(SymOper.U)
Uinv = nan;
elseif SymOper.conjugate
Uinv = conj(inv(SymOper.U));
else
Uinv = inv(SymOper.U);
end
Rinv = inv(SymOper.R);
tinv = SymOper.t*Rinv.';
SymOper = Oper(Rinv, Uinv,tinv,'conjugate',SymOper.conjugate,...
'antisymmetry',SymOper.antisymmetry,...
'strict_eq',SymOper.strict_eq...
);
end
