function [basic_eq,U_eq]=eq(SymOper1,SymOper2)
m = length(SymOper1);
n = length(SymOper2);
if m == 0 && n == 0
basic_eq =true;
U_eq = true;
elseif m == 0 || n == 0
basic_eq = false;
U_eq = false;
elseif m == 1 && n ==1
R_eq = Oper.allclose(SymOper1.R, SymOper2.R);
basic_eq = R_eq && SymOper1.conjugate ==SymOper2.conjugate && SymOper1.antisymmetry == SymOper2.antisymmetry && isequal(sym(SymOper1.t),  sym(SymOper2.t));
if basic_eq && (SymOper1.strict_eq || SymOper2.strict_eq)
if isnan(SymOper1.U) && isnan(SymOper2.U)
U_eq = true;
elseif xor(isnan(SymOper1.U) , isnan(SymOper2.U))
U_eq = false;
else
SymOper3 = (SymOper1.inv()*SymOper2);
[prop, coeff] = Oper.prop_to_id(SymOper3.U);
U_eq = (prop && Oper.isclose(abs(coeff), 1));
end
else
U_eq = true;
end
elseif m == 1
basic_eq(n) = false;
U_eq(n) = false;
for i = 1:n
[basic_eq(i),U_eq(i)] = eq(SymOper1,SymOper2(i));
end
elseif n == 1
basic_eq(m) = false;
U_eq(m) = false;
for i = 1:m
[basic_eq(i),U_eq(i)] = eq(SymOper2,SymOper1(i));
end
elseif m == n
basic_eq = false(size(SymOper1));
U_eq = false(size(SymOper1));
SymOper1 = sort(SymOper1);
SymOper2 = sort(SymOper2);
for i =1:m
[basic_eq(i),U_eq(i)] = eq(SymOper2(i),SymOper1(i));
end
else
basic_eq = false;
U_eq = false;
end
end
