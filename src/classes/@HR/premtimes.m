function H_hr = premtimes(A,B)
% PREMTIMES Left matrix multiplication for HR objects
%
%   H_HR = PREMTIMES(A,B) implements left multiplication (B*A) for HR objects
%
%   Inputs:
%       A - HR object
%       B - Numeric or symbolic matrix
%   Output:
%       H_hr - Resulting HR object
%
%   Notes:
%       - Only supports case where A is HR and B is numeric/symbolic
%       - Multiplies B from left to all hopping matrices
%       - Preserves HR object structure
if isa(A,'HR') && ~isa(B,'HR')
H_hr1 = A;
H_hr = H_hr1;
if isa(B,'sym')
for i = 1:H_hr.NRPTS
H_hr.HcoeL(:,:,i) = B*H_hr.HcoeL(:,:,i);
end
elseif isa(B,'numeric')
for i = 1:H_hr.NRPTS
H_hr.HnumL(:,:,i) = B*H_hr.HnumL(:,:,i) ;
end
else
end
end
end
