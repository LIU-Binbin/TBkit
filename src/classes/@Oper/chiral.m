function SymOper = chiral(realspace_dim, U,propArgs)
arguments
realspace_dim  double {mustBeInteger} =3;
U double =nan;
propArgs.?Oper;
propArgs.antisymmetry = true;
end
propertyCell = namedargs2cell(propArgs);
U_c = U;
R_c = eye(realspace_dim);
SymOper = Oper(R_c,U_c,propertyCell{:});
end
