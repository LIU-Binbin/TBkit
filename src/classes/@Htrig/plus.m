function C = plus(A,B)
if isa(A,'Htrig') && isa(B,'Htrig')
C = A;
mat_cell = mat2cell(reshape(B.HcoeL,B.Basis_num,B.Basis_num*B.Kinds).',repmat(B.Basis_num,[1 B.Kinds]));
Var_cell = mat2cell(sym(ones(1,B.Kinds)),1);
k_cell = mat2cell(B.HsymL_trig,1);
C = C.setup(Var_cell,k_cell,mat_cell);
elseif isa(A,'Htrig') && ~isa(B,'Htrig')
if isa(B,'Term')
C = A;
elseif isa(B,'Trig')
C = A;
C.Trig_list = C.Trig_list+B;
for i =1:length(B)
C = C.setup_rough(B(i).symbolic_polynomial,B(i).pauli_mat);
end
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
SymPoly = B(i,j);
A = A.setup_rough(SymPoly,tempmat);
end
end
end
C = A;
else
C = A;
end
elseif ~isa(A,'Htrig') && isa(B,'Htrig')
if isa(A,'Term')
C = B;
elseif isa(A,'Trig')
C = B;
C.Trig_list = C.Trig_list+A;
for i = 1:length(A)
C = C.setup_rough(A(i).symbolic_polynomial,A(i).pauli_mat);
end
else
C = B;
end
else
C = 0;
end
end
