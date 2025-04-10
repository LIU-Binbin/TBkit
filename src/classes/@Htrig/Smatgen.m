function [H_htrig,Smat] = Smatgen(H_htrig,R,Accuracy)
%SMATGEN  Generate symbolic transformation matrix for trig basis.
%
%   [H_htrig, Smat] = SMATGEN(H_htrig, R)
%   constructs the symbolic transformation matrix Smat for the trig
%   functions in H_htrig after a coordinate transformation defined by R.
%
%   [H_htrig, Smat] = SMATGEN(H_htrig, R, Accuracy)
%   performs the same procedure but allows specification of numerical 
%   rounding accuracy when R is numeric. Default Accuracy is 6 (i.e., 1e-6).
%
%   This function is useful in the context of tight-binding Hamiltonians
%   represented symbolically using trigonometric basis functions. The
%   transformation matrix Smat can be used to project or remap the 
%   Hamiltonian to a rotated coordinate frame.
%
%   Inputs:
%       H_htrig   - An instance of the Htrig class, containing symbolic
%                   Hamiltonian trigonometric terms (HsymL_trig).
%       R         - Coordinate transformation matrix (1x1, 2x2, or 3x3).
%                   Can be symbolic or numeric. Defines how momentum
%                   components [k_x, k_y, k_z] are linearly transformed.
%       Accuracy  - (Optional) Accuracy for rounding numerical Smat values
%                   if R is numeric. Default: 6 (i.e., round to 1e-6).
%
%   Outputs:
%       H_htrig   - Updated Htrig object. Additional basis functions may be
%                   added to HsymL_trig if new ones arise from transformation.
%       Smat      - Transformation matrix (symbolic or numeric), such that:
%                       H_transformed(k') = Smat * original_basis
%
%   Behavior:
%       - Substitutes new k vector (k') into the Hamiltonian.
%       - Extracts coefficients in terms of existing or extended basis.
%       - Builds a transformation matrix Smat mapping new H(k') to old basis.
%       - If R is numeric, result is rounded to given Accuracy using roundn.
%
%   Example:
%       syms kx ky kz real
%       H_htrig = Htrig();        % Create symbolic Htrig object
%       R = [0 -1; 1 0];          % 90-degree rotation in k-space
%       [H_htrig, Smat] = Smatgen(H_htrig, R, 6);
%
%   See also: Htrig, coeffs, simplify, subs
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
