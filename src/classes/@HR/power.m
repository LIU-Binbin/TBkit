function H_hr = power(H_hr,b)
% POWER Element-wise power operation for HR object
%
%   H_HR = POWER(H_HR,B) raises HR matrix elements to power B
%
%   Inputs:
%       H_hr - HR object to modify
%       b - Power to raise elements
%   Output:
%       H_hr - Modified HR object
%
%   Notes:
%       - Operates on both numeric and symbolic coefficients
%       - Performs element-wise power (.^) operation
%       - Processes all hopping terms
for i = 1:H_hr.NRPTS
H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i).^b;
H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i).^b;
end
end
