function SymOper = mirror(axis, U, spin,options,propArgs)
arguments
    axis  = nan;
    U  =nan;
    spin double = nan;
    options.sym = false;
    propArgs.?Oper;
end
propertyCell = namedargs2cell(propArgs);
U_mirror = U;
if ~all(all(isnan(U))) && ~all(all(all(isnan(spin))))
    raise ValueError('Only one of `U` and `spin` may be provided.');
end
axis = axis/norm(axis);
axis = reshape(axis,length(axis),1);
R_mirror = eye(length(axis)) - 2 *(axis*axis.');
if ~isnan(spin)
    if length(axis) == 2
        axis(3) = 0;
    end
    U_mirror = Oper.spin_rotation(pi * axis, spin);
end
if isa(R_mirror,'sym')
    R_mirror = simplify(R_mirror);
end
if isa(U_mirror,'sym')
    U_mirror = simplify(U_mirror);
end
SymOper = Oper(R_mirror,U_mirror,propertyCell{:});
end
