function C = isclose(A,B)
if isa(A,'sym') ||  isa(B,'sym')
    C = logical(simplify(sym(A))==simplify(sym(B)));
else
    C = (A-B).^2<1e-12 ;
end
end
