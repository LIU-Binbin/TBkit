function C = mrdivide(A,B)
A = contract(A);
B = contract(B);
if A == B
C = A(1).coe/B(1).coe;
else
C = 0;
end
end
