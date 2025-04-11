function H_hr = uminus(H_hr)
%UMINUS Unary minus for HR objects
%
%   H_HR = UMINUS(H_HR) Negates all Hamiltonian matrix elements
%
%   Input:
%       H_hr - HR object to negate
%
%   Output:
%       H_hr - HR object with negated elements
%
%   See also HR, UPLUS, MINUS
for i = 1:H_hr.NRPTS
H_hr.HcoeL(:,:,i) = -H_hr.HcoeL(:,:,i);
H_hr.HnumL(:,:,i) = -H_hr.HnumL(:,:,i);
end
end
