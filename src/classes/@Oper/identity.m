function SymOper = identity(dim, shape,propArgs)
arguments
    dim  double {mustBeInteger} =3;
    shape = nan;
    propArgs.?Oper;
end
propertyCell = namedargs2cell(propArgs);
if ~isnan(shape)
    U_ = eye(shape);
else
    U_ = nan;
end
SymOper = Oper(eye(dim),U_,propertyCell{:});
SymOper.isIdentity = true;

end
