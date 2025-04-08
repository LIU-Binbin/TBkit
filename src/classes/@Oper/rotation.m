function SymOper = rotation(angle, axis, inversion, U, spin,options,propArgs)
arguments
angle  ;
axis  = nan;
inversion logical = false;
U  = nan;
spin  = nan;
options.sym = false;
options.rightorleft {mustBeMember(options.rightorleft,{'left','right'})}= 'right';
propArgs.?Oper;
end
propertyCell = namedargs2cell(propArgs);
if isa(U,'double')
if ~all(all(isnan(U))) && ~all(all(all(isnan(spin))))
raise ValueError('Only one of `U` and `spin` may be provided.');
end
elseif isa(U,'pauli_matrix')
U = double(U);
end
if options.sym
if isa(angle,'sym')
else
angle = sym(2 * pi * angle);
end
U_PG = sym(U);
else
if isa(angle,'sym')
else
angle = 2 * pi * angle;
end
U_PG = U;
end
if strcmp(options.rightorleft,'left')
rightorleft = -1;
else
rightorleft = 1;
end
if isnan(axis)
R_PG = ([[cos(angle), sin(angle)];
[-sin(angle), cos(angle)]]);
if ~isnan(spin)
U_PG = Oper.spin_rotation(angle * ([0, 0, 1]), spin);
end
elseif length(axis) == 3
n = angle * axis / norm(axis);
R_PG = Oper.spin_rotation(n, rightorleft*Oper.L_matrices(3, 1));
if inversion
R_PG = -R_PG;
else
end
if ~isnan(spin)
if ~isvector(spin)
U_PG = Oper.spin_rotation(n, spin);
else
end
end
else
raise ValueError('`axis` needs to be `None` or a 3D vector.')
end
SymOper = Oper(R_PG,U_PG,propertyCell{:});
end
