function [prop,Coeff] = prop_to_id(A)
if Oper.isclose(A(1,1), 0)
if matallclose(A, zeros(length(A)))
prop = true;
Coeff = 0;
else
prop = false;
Coeff = 0;
end
else
if Oper.allclose( A /A(1,1), eye(length(A)))
prop = true;
Coeff = A(1,1);
else
prop = false;
Coeff = 0;
end
end
end
