function H_hr = transpose(H_hr)
%TRANSPOSE Overloaded transpose for HR objects
%
%   H_HR = TRANSPOSE(H_HR) computes the non-conjugate transpose
%   of the Hamiltonian matrices in an HR object.
%
%   Input:
%       H_hr - HR object to transpose
%
%   Output:
%       H_hr - Transposed HR object
%
%   See also HR, CTRANSPOSE, PAGETRANSPOSE
H_hr.HcoeL = pagetranspose(H_hr.HcoeL);
H_hr.HnumL = pagetranspose(H_hr.HnumL);
end
