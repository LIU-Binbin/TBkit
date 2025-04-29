function C = isclose(A,B)
%ISCLOSE Check if values are numerically close within tolerance
%
%   Syntax:
%       C = isclose(A,B)
%
%   Description:
%       Compares arrays A and B element-wise with a fixed tolerance (1e-12).
%       Returns logical array indicating where values are nearly equal.
%
%   Inputs:
%       A - First input array
%       B - Second input array
%
%   Output:
%       C - Logical array where |A-B|² < 1e-12
%
%   See also: allclose
C = (A-B).^2<1e-12 ;
end