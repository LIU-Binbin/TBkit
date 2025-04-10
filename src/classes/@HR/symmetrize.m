function H_hr = symmetrize(H_hr,SymOper,options)
% SYMMETRIZE Apply symmetry operations to HR object
%
%   H_HR = SYMMETRIZE(H_HR,SYMOPER,OPTIONS) symmetrizes Hamiltonian
%
%   Inputs:
%       H_hr - HR object
%       SymOper - Symmetry operator [default: Oper()]
%       options.generator - Generator mode [default: false]
%   Output:
%       H_hr - Symmetrized HR object
%
%   Notes:
%       - Under development
%       - Will handle space group symmetries
%       - Supports generator mode
arguments
H_hr HR;
SymOper Oper = Oper();
options.generator = false;
end
end
