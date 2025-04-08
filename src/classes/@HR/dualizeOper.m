function [H_hr,VectorDistMat] = dualizeOper(H_hr,SymOper)
Accuracy = 1e-6;
NRPTS_ = H_hr.NRPTS;
Rf = double(Oper.Rc2Rf(inv(SymOper.R),H_hr.Rm));
U = SymOper.U;
invRf = inv(double(Rf));
invU = inv(U);
[ml_cell,ij_list] = Tij2lm(H_hr,Rf);
vectorList = double(H_hr.vectorL(:,1:H_hr.Dim));
ij_list_full = double(H_hr.vectorL(:,H_hr.Dim+1:H_hr.Dim+2));
VectorDistMat = zeros(NRPTS_,NRPTS_);
for in = 1:NRPTS_
ij = ij_list_full(in,:);
Rvector = vectorList(in,:);
[~,k] = ismember(ij,ij_list,'rows');
ml_list = ml_cell{k};
i = ij(1);
j = ij(2);
for kn  = 1 : size(ml_list,1)
T_ij__ml = ml_list(kn,3:5);
m =  ml_list(kn,1);
l =  ml_list(kn,2);
vector_tmp_oppo = (Rvector - T_ij__ml)*invRf;
tmp_vector = ([vector_tmp_oppo,[l,m]]);
[~,jn] = ismember(tmp_vector,H_hr.vectorL,'rows');
if jn == 0
H_hr = H_hr.add_empty_one(tmp_vector);
jn = H_hr.NRPTS;
U_tmp = U(j,m)*invU(l,i);
U_tmp_r = real(U_tmp);
U_tmp_i = imag(U_tmp);
if abs(U_tmp_r) < Accuracy
U_tmp_r = 0 ;
end
if abs(U_tmp_i) < Accuracy
U_tmp_i = 0 ;
end
VectorDistMat(jn,in) = U_tmp_r + 1i * U_tmp_i ;
end
U_tmp = U(i,l)*invU(m,j);
U_tmp_r = real(U_tmp);
U_tmp_i = imag(U_tmp);
if abs(double(U_tmp_r)) < Accuracy
U_tmp_r = 0 ;
end
if abs(double(U_tmp_i)) < Accuracy
U_tmp_i = 0 ;
end
VectorDistMat(in,jn) = U_tmp_r + 1i * U_tmp_i ;
end
end
end
