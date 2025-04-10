function H_hr_bk = subsOper(H_hr,SymOper)
% SUBSOPER Apply symmetry operator substitution
%
%   H_HR_BK = SUBSOPER(H_HR,SYMOPER) applies symmetry operation
%
%   Inputs:
%       H_hr - Htrig object
%       SymOper - Symmetry operator (Oper object)
%   Output:
%       H_hr_bk - Transformed HR object
%
%   Notes:
%       - Handles unitary transformations
%       - Supports conjugation and antisymmetry
%       - Initializes object if needed
%       - Returns backup of original
arguments
H_hr Htrig;
SymOper Oper = Oper();
end
if isequal(sym(zeros(size(H_hr.HcoeL))),H_hr.HcoeL)
H_hr.coe = false;
else
H_hr.coe = true;
end
if ~H_hr.coe
H_hr = H_hr.init();
H_hr = H_hr.hermitize();
end
if length(SymOper) == 1
if ~SymOper.conjugate && ~SymOper.antisymmetry && isequal(SymOper.R,eye(3))
return;
end
[H_hr_bk,H_hr]  = H_hr.applyR(inv(SymOper.R));
if isnan(SymOper.U)
end
H_hr_bk  = H_hr_bk.applyU(SymOper.U,SymOper.conjugate,SymOper.antisymmetry);
end
end
