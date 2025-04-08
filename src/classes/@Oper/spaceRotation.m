function SymOper = spaceRotation(angle, axis,t, inversion, U, spin,options,propArgs)
arguments
angle double ;
axis double = nan;
t double =[0,0,0];
inversion logical = false;
U double =nan;
spin double = nan;
options.sym = false;
propArgs.?Oper;
end
propertyCell = namedargs2cell(propArgs);
U_SG = U;
t_SG = t;
if ~all(all(isnan(U))) && ~all(all(all(isnan(spin))))
raise ValueError('Only one of `U` and `spin` may be provided.');
end
angle = 2 * pi * angle;
if isnan(axis)
R_SG = ([[cos(angle), -sin(angle)];
[sin(angle), cos(angle)]]);
if ~isnan(spin)
U_SG = Oper.spin_rotation(angle * ([0, 0, 1]), spin);
end
elseif len(axis) == 3
n = angle * axis / norm(axis);
R_SG = Oper.spin_rotation(n, L_matrices(d=3, l=1));
if inversion
R_SG = -R_SG;
else
end
if ~isnan(spin)
U_SG = spin_rotation(n, spin);
end
else
raise ValueError('`axis` needs to be `None` or a 3D vector.')
end
SymOper = Oper(R_SG,U_SG,t_SG,propertyCell{:});
end
