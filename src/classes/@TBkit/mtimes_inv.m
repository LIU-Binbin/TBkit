function C = mtimes_inv(B,A)
%MTIMES_INV Reverse matrix multiplication for HK objects
%   C = MTIMES_INV(B, A) performs reverse matrix multiplication B*A
%   where B is HK object and A is numeric/symbolic.
%
%   Inputs:
%       B - HK object
%       A - Numeric or symbolic matrix
%
%   Output:
%       C - Result of multiplication operation
%
%   See also HK, MTIMES
if ~isa(A,'HK') && isa(B,'HK')
    C = B;
    if isa(A,'double')
        for i = 1:C.Kinds
            C.HcoeL(:,:,i) = A*C.HcoeL(:,:,i);
            C.HnumL(:,:,i) = A*C.HnumL(:,:,i);
        end
    elseif isa(A,'sym')
        for i = 1:C.Kinds
            C.HcoeL(:,:,i) = A*C.HcoeL(:,:,i);
        end
    else
    end
end
end
