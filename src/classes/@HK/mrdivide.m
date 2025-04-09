function C = mrdivide(A,B)
%MRDIVIDE Overloaded matrix right division for HK objects
%
% Syntax:
%   C = A / B       % Operator form
%   C = mrdivide(A,B) % Functional form
%
% Inputs:
%   A - HK object
%   B - Scalar divisor
%
% Output:
%   C - Scaled HK object
%
% Description:
%   Divides all Hamiltonian coefficients by scalar B:
%   - Numeric coefficients (HnumL) if present
%   - Symbolic coefficients (HcoeL)
%
% Note:
%   Only implemented for HK/scalar case
%   Preserves term structure and other properties
%
% Example:
%   Hk_scaled = Hk / 2; % Halves all coefficients
if isa(A,'HK') && isa(B,'HK')
elseif isa(A,'HK') && ~isa(B,'HK')
    C = A;
    C.HcoeL =  C.HcoeL/B;
else
end
end
