function C = mtimes(A,B)
%MTIMES Overloaded matrix multiplication for HK objects
%
% Syntax:
%   C = A * B       % Operator form
%   C = mtimes(A,B) % Functional form
%
% Inputs:
%   A - HK object or numeric/symbolic matrix
%   B - HK object or numeric/symbolic matrix
%
% Output:
%   C - Result of multiplication
%
% Description:
%   Handles three cases of HK multiplication:
%   1. Scalar * HK: Scales all coefficients
%   2. HK * Scalar: Scales all coefficients
%   3. HK * Matrix: Right-multiplies each coefficient matrix
%   4. Matrix * HK: Left-multiplies each coefficient matrix
%
% Note:
%   HK*HK multiplication not currently implemented
%   Preserves term structure when scaling
%
% Example:
%   Hk_scaled = 3 * Hk; % Triples all coefficients
%   Hk_transformed = U * Hk; % Transforms basis
if isa(A,'HK') && isa(B,'HK')
    C = A;
elseif isa(A,'HK') && ~isa(B,'HK')
    C =A;
    if isa(B,'double')
        for i = 1:C.Kinds
            C.HcoeL(:,:,i) = C.HcoeL(:,:,i)*B;
            C.HnumL(:,:,i) = C.HnumL(:,:,i)*B;
        end
    elseif isa(B,'sym')
        for i = 1:C.Kinds
            C.HcoeL(:,:,i) = C.HcoeL(:,:,i)*B;
        end
    else
    end
elseif ~isa(A,'HK') && isa(B,'HK')
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
else
    C = 0;
end
end
