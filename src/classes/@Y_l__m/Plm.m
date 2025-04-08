function PlmExpr = Plm(l,m,seed,options)
arguments
l double{mustBeNonnegative,mustBeInteger};
m double{mustBeInteger}
seed sym{} = (sym('theta','real'));
options.ClosedForm = false;
options.triangle = true;
end
optionsCell = namedargs2cell(options);
if abs(m)>l
ME = MException('TBkit:Y_l__m:WrongInput', ...
'Variable m:
throw(ME);
end
if strcmp(string(seed),"theta") && ~options.triangle
x = sym('x','real');
y = (1-x^2)^0.5;
elseif options.triangle
x  = cos(seed);
y  = sin(seed);
else
x = seed;
y = (1-x^2)^0.5;
end
if options.ClosedForm
if  options.triangle
PlmExpr_pre = (-1)^m * 2^l * (y)^m;
else
PlmExpr_pre = (-1)^m * 2^l * ((1-x^2)^(1/2))^m;
end
binomialfunc = @nchoosek;
k = m;
PlmExpr_tail = factorial(k)/factorial(k-m)*x^(k-m)*binomialfunc(sym(l),k)*binomialfunc((l+k-1)/2,sym(l));
for k = m+1:l
PlmExpr_tail =PlmExpr_tail + factorial(k)/factorial(k-m)*x^(k-m)*binomialfunc(sym(l),k)*binomialfunc((l+k-1)/2,sym(l));
end
PlmExpr = simplify(PlmExpr_pre*simplify(PlmExpr_tail));
return;
end
if m < 0
PlmExpr = (-1)^(-m) * factorial(l+m)/factorial(l-m) * Y_l__m.Plm(-m,l,seed,optionsCell{:});
return;
end
if l==0 && m ==0
PlmExpr = 1;
return;
elseif l == m
PlmExpr = (-1)^l*Y_l__m.DFactorial(2*l-1)*y^l;
return;
elseif m == l-1
PlmExpr = (2*m+1)*x*Y_l__m.Plm(m,m,seed,optionsCell{:});
else
PlmExpr = (2*l-1)/(l-m) *x * Y_l__m.Plm(l-1,m,seed,optionsCell{:})...
+ (1-l-m)/(l-m) * Y_l__m.Plm(l-2,m,seed,optionsCell{:});
end
end
