function Factorlist_parity = factorlist_parity(Degree)
%FACTORLIST_PARITY Generate parity factors for k·p terms
%
% Syntax:
%   Factorlist_parity = factorlist_parity(Degree)
%
% Input:
%   Degree - Maximum polynomial degree
%
% Output:
%   Factorlist_parity - Array of parity factors (+1/-1) for each term
%
% Description:
%   Generates parity factors (+1 for even, -1 for odd) for all terms in
%   a k·p expansion up to specified degree. Used for time-reversal symmetric
%   operations and other transformations requiring parity consideration.
%
% Algorithm:
%   Terms are classified by polynomial degree:
%   - Even degree terms get +1
%   - Odd degree terms get -1
%   Uses HK.orderlist to maintain consistent term ordering
%
% Example:
%   factors = factorlist_parity(2) % Returns [1 -1 -1 -1 1 1 1 1 1 1]
Nterm = nchoosek(4+Degree-1,Degree);
Factorlist_parity = zeros(1,Nterm);
for i =0:Degree
    if mod(i, 2) == 0
        Factorlist_parity(HK.orderlist(i)) = 1;
    else
        Factorlist_parity(HK.orderlist(i)) = -1;
    end
end
end
