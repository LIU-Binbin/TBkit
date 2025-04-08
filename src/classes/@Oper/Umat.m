function U = Umat(SymOper1,SymOper2,options)
arguments
SymOper1
SymOper2
options.disp = true;
end
if isa(SymOper1,'Oper')
A =SymOper1.U;
else
A =SymOper1;
end
if isa(SymOper2,'Oper')
B =SymOper2.U;
else
B =SymOper2;
end
[QA,DA]=eig(A);
[QB,DB]=eig(B);
k = DA/DB;
U = sym((QB/QA)/sqrtm(k));
if options.disp
disp(('U * A * U^* = B'));
sym(A)
sym(B)
sym(U)
end
end
