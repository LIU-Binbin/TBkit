function H_htrig = rotation(H_htrig,Rotation)
if nargin <2
Rotation = inv(sym(H_htrig.Rm));
end
if ischar(Rotation)
if strcmp(Rotation,'auto')
Rotation = inv(sym(H_htrig.Rm));
end
end
VarUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
k = Rotation * VarUsing.';
for i = 1:H_htrig.Dim
VarUsingCell{i} = VarUsing(i);
kCell{i} = k(i);
end
H_htrig = H_htrig.subs(VarUsingCell,kCell);
end
