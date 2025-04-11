function [H_hr,Smat] = Smatgen(H_hr,R,Accuracy)
% SMATGEN Generate S-matrix for HR object
%
%   [H_HR,Smat] = SMATGEN(H_HR,R,ACCURACY) generates S-matrix transformation
%
%   Inputs:
%       H_hr - Htrig object
%       R - Transformation matrix
%       Accuracy - Rounding accuracy [default: 6]
%   Outputs:
%       H_hr - Modified HR object
%       Smat - Generated S-matrix
%
%   Notes:
%       - Handles symbolic transformations
%       - Supports 1D, 2D and 3D systems
%       - Maintains basis consistency
%       - Can switch between numeric/symbolic modes
arguments
H_hr Htrig;
R  ;
Accuracy double = 6;
end
if isa(R,'sym')
H_hr.coe = true;
else
H_hr.coe = false;
end
BASIS_NUM = H_hr.Basis_num;
HsymC = H_hr.HsymL_trig;
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
[coeff_trig,~] = coeffs(HsymC(k),H_hr.HsymL_trig_bk);
for i = 1:numel(coeff_trig)
tmp_label = contains(string(coeff_trig),H_hr.seeds);
if sum(tmp_label)
[~,coeff_trig_list] = coeffs(coeff_trig(i));
for j = 1:numel(coeff_trig_list)
tmp_label2 = contains(string(coeff_trig_list(j)),H_hr.seeds);
if sum(tmp_label2)
H_hr = H_hr.find_HsymL_trig_bk(coeff_trig_list(j));
end
end
end
end
end
for i = 1:numel(HsymC)
[~,B] = coeffs(HsymC(i),H_hr.HsymL_trig_bk);
for k = 1:length(B)
Kind = H_hr.k_symbol2Kind(B(k));
if isempty(Kind)
Kind = H_hr.Kinds+1;
H_hr.HsymL_trig(Kind) = B(k);
H_hr.HcoeL(:,:,Kind) = sym(zeros(BASIS_NUM,BASIS_NUM,1));
H_hr.HnumL(:,:,Kind)  = (zeros(BASIS_NUM,BASIS_NUM,1));
end
end
end
HsymC_bk = H_hr.HsymL_trig;
Smat =sym(zeros(numel(HsymC_bk)));
for i = 1:numel(HsymC)
[A,B] = coeffs(HsymC(i),H_hr.HsymL_trig_bk);
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
if ~H_hr.coe
Smat = roundn(double(Smat),-Accuracy);
end
end
