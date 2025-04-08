function SymOper = E(SymOper)
dim = length(SymOper.R);
R_E = eye(dim) ;
if isnan(SymOper.U)
U_E = SymOper.U;
else
U_E = eye(length(SymOper.U));
end
SymOper = Oper(R_E, U_E,SymOper.t,'conjugate',SymOper.conjugate,...
'antisymmetry',SymOper.antisymmetry,...
'strict_eq',SymOper.strict_eq...
);
end
