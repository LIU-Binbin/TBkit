function C = allclose(A,B)
if all(all(Oper.isclose(A,B)))
C = true;
else
C = false;
end
end
