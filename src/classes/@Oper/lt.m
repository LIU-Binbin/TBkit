function result = lt(SymOper1,SymOper2)
Rs = SymOper1.R;
Ro = SymOper2.R;
identity = eye(length(Rs));
if SymOper1.conjugate ~= SymOper2.conjugate || ...
SymOper1.antisymmetry ~= SymOper2.antisymmetry
if SymOper1.conjugate < SymOper2.conjugate
result =true;
elseif SymOper1.conjugate == SymOper2.conjugate
if SymOper1.antisymmetry < SymOper2.antisymmetry
result =true;
else
result =false;
end
else
result =false;
end
elseif xor(Oper.allclose(Rs,identity), Oper.allclose(Ro,identity))
result = Oper.allclose(Rs,identity);
else
L = Rs(:) < Ro(:);
B = ~Oper.isclose(Rs(:),Ro(:));
for i =1:length(B)
if B(i)
result = L(i);
return;
end
end
result  = false;
end
end
