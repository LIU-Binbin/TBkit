function SymOper = inversion(realspace_dim, U,quantumL,propArgs)
arguments
realspace_dim  double {mustBeInteger} =3;
U double =nan;
quantumL = nan;
propArgs.?Oper;
end
propertyCell = namedargs2cell(propArgs);
U_inv = U;
if ~isnan(quantumL)
end
R_inv = -eye(realspace_dim);
SymOper = Oper(R_inv,U_inv,propertyCell{:});
end
