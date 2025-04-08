function C = plus(A,B)
if isa(A,'HK') && isa(B,'HK')
H_hk1 = A;
H_hk2 = B;
if H_hk1.Basis_num ~= H_hk2.Basis_num
error('basis num differ');
end
if H_hk1.Degree ~= H_hk2.Degree
error('Degree differ');
end
C = H_hk1;
C.HcoeL = C.HcoeL + H_hk2.HcoeL;
elseif isa(A,'HK') && ~isa(B,'HK')
if isa(B,'Term')
C = A;
if contains(C.Type,'kp')
C.Term_to_save = C.Term_to_save + B;
else
C.Term_to_save = B;
end
for i = 1:length(B)
C = C.setup_rough(B(i).symbolic_polynomial,B(i).pauli_mat);
end
elseif isa(B,'Trig')
C = A;
C.Trig_to_save = C.Trig_to_save+B.symbolic_polynomial*double(B.pauli_mat);
elseif isa(B,'sym')
basis_num = A.Basis_num;
if basis_num~= length(B) || size(B,2)  ~=  size(B,2)
error('Basis wrong!')
end
for i = 1:basis_num
for j = 1: basis_num
if B(i,j)~=sym(0)
tempmat = zeros( basis_num);
tempmat(i,j) =1 ;
A = A.setup_rough(B(i,j),tempmat);
end
end
end
C =A;
else
C = A;
end
elseif ~isa(A,'HK') && isa(B,'HK')
if isa(A,'Term')
C = B;
if contains(C.Type,'kp')
C.Term_to_save = C.Term_to_save + A;
else
C.Term_to_save = A;
end
for i = 1:length(A)
C = C.setup_rough(A(i).symbolic_polynomial,A(i).pauli_mat);
end
elseif isa(A,'Trig')
C = B;
C.Trig_to_save = C.Trig_to_save + A.symbolic_polynomial*double(A.pauli_mat);
else
C = B;
end
else
C = 0;
end
end
