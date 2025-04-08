function SymOper = particle_hole(realspace_dim, U,propArgs)
arguments
realspace_dim  double {mustBeInteger} =3;
U double =nan;
propArgs.?Oper;
propArgs.conjugate = true;
propArgs.antisymmetry = true;
end
propertyCell = namedargs2cell(propArgs);
U_ph = U;
R_ph = eye(realspace_dim);
SymOper = Oper(R_ph,U_ph,propertyCell{:});
end
