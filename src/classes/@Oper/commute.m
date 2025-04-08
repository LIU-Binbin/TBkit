function [TrueOrNot,result] = commute(SymOper1,SymOper2)
U1 = SymOper1.U;
U2 = SymOper2.U;
result = U1*U2-U2*U1;
if isa(U1,'sym')||isa(U2,'sym')
zeromat = sym(zeros(size(U1)));
else
zeromat = zeros(size(U1));
end
if isequal(result,zeromat)
TrueOrNot = true;
else
TrueOrNot = false;
end
end
