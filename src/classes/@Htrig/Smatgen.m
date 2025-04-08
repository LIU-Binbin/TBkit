function [H_htrig,Smat] = Smatgen(H_htrig,R,Accuracy)
arguments
H_htrig Htrig;
R  ;
Accuracy double = 6;
end
if isa(R,'sym')
coe_label = true;
else
coe_label = false;
end
BASIS_NUM = H_htrig.Basis_num;
HsymC = H_htrig.HsymL_trig;
syms k_x k_y k_z real;
if length(R) == 3
k = R*[k_x; k_y; k_z];
HsymC = subs(HsymC,[k_x k_y k_z],k.');
elseif length(R) == 2
k = R*[k_x; k_y];
HsymC = subs(HsymC,[k_x k_y],k.');
elseif length(R) == 1
k = R*[k_x];
HsymC = subs(HsymC,[k_x],k.');
end
HsymC = expand(simplify(HsymC));
for k = 1:length(HsymC)
[coeff_trig,~] = coeffs(HsymC(k),H_htrig.HsymL_trig_bk);
for i = 1:numel(coeff_trig)
tmp_label = contains(string(coeff_trig),H_htrig.seeds);
if sum(tmp_label)
[~,coeff_trig_list] = coeffs(coeff_trig(i));
for j = 1:numel(coeff_trig_list)
tmp_label2 = contains(string(coeff_trig_list(j)),H_htrig.seeds);
if sum(tmp_label2)
H_htrig = H_htrig.find_HsymL_trig_bk(coeff_trig_list(j));
end
end
end
end
end
for i = 1:numel(HsymC)
[~,B] = coeffs(HsymC(i),H_htrig.HsymL_trig_bk);
for k = 1:length(B)
Kind = H_htrig.k_symbol2Kind(B(k));
if isempty(Kind)
Kind = H_htrig.Kinds+1;
H_htrig.HsymL_trig(Kind) = B(k);
H_htrig.HcoeL(:,:,Kind) = sym(zeros(BASIS_NUM,BASIS_NUM,1));
H_htrig.HnumL(:,:,Kind)  = (zeros(BASIS_NUM,BASIS_NUM,1));
end
end
end
HsymC_bk = H_htrig.HsymL_trig;
Smat =sym(zeros(numel(HsymC_bk)));
for i = 1:numel(HsymC)
[A,B] = coeffs(HsymC(i),H_htrig.HsymL_trig_bk);
for k = 1:length(B)
tempSym = B(k);
for l  = 1:numel(HsymC_bk)
if isequal(tempSym,HsymC_bk(l))
Smat(l,i)=A(k);
break;
end
end
end
end
if ~coe_label
Smat = roundn(double(Smat),-Accuracy);
end
end
