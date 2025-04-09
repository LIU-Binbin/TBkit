function C = mtimes_inv(B,A)
%MTIMES_INV Inverse multiplication ordering for HK objects
%
% Syntax:
%   C = mtimes_inv(B,A)
%
% Inputs:
%   B - HK object
%   A - Numeric/symbolic matrix
%
% Output:
%   C - Transformed HK object
%
% Description:
%   Specialized left-multiplication that:
%   1. Applies matrix A to each coefficient matrix
%   2. Handles both numeric (HnumL) and symbolic (HcoeL) coefficients
%   3. Preserves all other Hamiltonian properties
%
% Note:
%   Equivalent to A*B when B is HK object
%   More efficient than direct multiplication for some operations
%
% Example:
%   Hk_transformed = mtimes_inv(U,Hk); % Same as U*Hk
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
