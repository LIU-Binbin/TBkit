function H_htrig = setup_rough(H_htrig,symbolic_polynomial,pauli_mat,silence)
if nargin < 4
silence = true;
end
[coeff_trig,symvar_list_trig,H_htrig] = H_htrig.split_sym_eq((simplify(symbolic_polynomial,'IgnoreAnalyticConstraints',true)));
nc = length(coeff_trig);
if nc >0
switch H_htrig.Type
case 'mat'
for i =1:nc
H_htrig = H_htrig.set_hop(coeff_trig(i)*sum(double(pauli_mat),3),symvar_list_trig(i,:));
end
case 'list'
Mtmp  = sum(double(pauli_mat),3);
LabelnonZero = find(Mtmp~=0);
nLabelnonZero = length(LabelnonZero);
[iL,jL] = ind2sub(size(Mtmp),LabelnonZero);
Mtmp = Mtmp(LabelnonZero);
Mlist = Mtmp(:);
Mtmp_list = kron(Mlist,ones(nc,1));
iL = kron(iL,ones(nc,1,'sym'));
jL = kron(jL,ones(nc,1,'sym'));
symvar_list_trig = repmat(symvar_list_trig,[nLabelnonZero,1]);
coeff_trig = repmat(coeff_trig(:),[nLabelnonZero,1]);
H_htrig = H_htrig.set_hop(coeff_trig.*Mtmp_list,[symvar_list_trig,iL,jL]);
otherwise
for i =1:nc
k_cell{i} = symvar_list_trig(i);
mat_cell{i} = sum(double(pauli_mat),3);
Var_cell{i} = coeff_trig(i);
end
H_htrig = setup(H_htrig,Var_cell,k_cell,mat_cell,silence);
end
end
end
