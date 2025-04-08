function SymOper = time_reversal(realspace_dim, U, spin,propArgs)
arguments
realspace_dim  double {mustBeInteger} =3;
U double =nan;
spin double = nan;
propArgs.?Oper;
propArgs.conjugate = true;
end
propertyCell = namedargs2cell(propArgs);
if ~all(all(isnan(U))) && ~all(all(all(isnan(spin))))
raise ValueError('Only one of `U` and `spin` may be provided.');
end
if ~isnan(spin)
U_tr = Oper.spin_rotation(pi*[0, 1, 0], spin);
else
U_tr = U;
end
R_tr = eye(realspace_dim);
SymOper = Oper(R_tr,U_tr,propertyCell{:});
end
