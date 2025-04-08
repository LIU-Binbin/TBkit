function P = similarity(SymOper1,SymOper2,options)
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
[V,~]=jordan(A);
[U,~]=jordan(B);
P = (V/(U));
if options.disp
disp(('P * A * P^-1 = B'));
A
B
P
end
end
