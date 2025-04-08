function C = horzcat(A,B)
if isa(A,'Htrig') && isa(B,'Htrig')
H_htrig1 = A;
H_htrig2 = B;
H_htrig = A;
H_htrig.Basis = [H_htrig1.Basis;H_htrig2.Basis];
H_htrig.Basis_num = H_htrig1.Basis_num+H_htrig2.Basis_num;
H_htrig.HnumL = zeros(H_htrig.Basis_num,H_htrig.Basis_num,H_htrig.Kinds);
H_htrig.HcoeL = sym(H_htrig.HnumL);
for i = 1:H_htrig.Kinds
H_htrig.HnumL(:,:,i) = blkdiag(H_htrig1.HnumL(:,:,i) ,...
H_htrig2.HnumL(:,:,i));
H_htrig.HcoeL(:,:,i) = blkdiag(H_htrig1.HcoeL(:,:,i) ,...
H_htrig2.HcoeL(:,:,i));
end
H_htrig.Trig_to_save =sym(zeros(H_htrig.Basis_num,H_htrig.Basis_num));
else
end
C = H_htrig;
end
