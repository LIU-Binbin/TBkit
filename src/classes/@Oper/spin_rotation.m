function U = spin_rotation(n, s, roundint)
arguments
n ;
s double  =3/2;
roundint logical = false;
end
if length(s) == 1
J_mat = Oper.spin_matrices(s);
else
J_mat = s;
end
expmat = Oper.tensordot_naive(n,J_mat);
U = expm(1i*expmat);
if roundint
Ur = round(real(U));
if Oper.allclose(U,Ur)
U = round(Ur);
else
error('?');
end
end
end
