function YlmExpr = InnerY_l__m(l,m,seed1,seed2,options)
arguments
l double{mustBeNonnegative,mustBeInteger};
m double{mustBeInteger}
seed1 sym{} = (sym('theta','real'));
seed2 sym{} = (sym('phi','real'));
options.ClosedForm = false;
options.triangle = true;
end
optionsCell = namedargs2cell(options);
if options.triangle
phi  = seed2;
YlmExpr = (-1)^m/sqrt(2*sym(pi))*...
Y_l__m.Nlm(l,m,seed1,optionsCell{:})*exp(1i*m*phi);
return;
else
if strcmp(string(seed2),"phi") &&  strcmp(string(seed1),"theta")
x = sym('x','real');
y = sym('y','real');
z = sym('z','real');
r = sym('r','real');
else
x = seed1(1); y = seed1(2); z =seed1(3);
r = seed2;
end
YlmExpr = sqrt((2*l+1)*factorial(l-m)/factorial(l+m)/(4*sym(pi)))* ...
Y_l__m.Plm(l,m,z/r,optionsCell{:})*(1-(z/r)^2)^(-m/2)*(x/r+1i*y/r)^m;
YlmExpr = simplify(YlmExpr);
end
end
