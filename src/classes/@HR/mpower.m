function H_hr = mpower(A,B)
% MPOWER Overloaded power operator for HR objects
%
%   H_HR = MPOWER(A,B) implements matrix power operation for HR objects
%
%   Inputs:
%       A - HR object or numeric scalar
%       B - HR object or numeric scalar
%   Output:
%       H_hr - Resulting HR object after power operation
%
%   Notes:
%       - Only supports scalar numeric exponents
%       - Implements both A^B and B^A cases for HR objects
%       - Throws error for invalid input combinations

if isa(A,'HR') && isa(B,'numeric')
    if length(B) >1
        error('only support a number');
    end
    H_hr = A;
    c = B;
    for i = 1:H_hr.NRPTS
        if H_hr.coe
            H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i)^c;
        end
        if H_hr.num
            H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i)^c;
        end
    end
elseif isa(B,'HR') && isa(A,'numeric')
    if length(A) >1
        error('only support a number');
    end
    H_hr = B;
    c = A;
    for i = 1:H_hr.NRPTS
        if H_hr.coe
            H_hr.HcoeL(:,:,i) = c^H_hr.HcoeL(:,:,i);
        end
        if H_hr.num
            H_hr.HnumL(:,:,i) = c^H_hr.HnumL(:,:,i);
        end
    end
else
    error('wrong input');
end
end
