function Orderlist = orderlist(order)
%ORDERLIST Generate term indices for specific polynomial order
%
% Syntax:
%   Orderlist = orderlist(order)
%
% Input:
%   order - Polynomial order (0 for constants)
%
% Output:
%   Orderlist - Vector of term indices for the specified order
%
% Description:
%   Computes the range of term indices corresponding to a given polynomial
%   order in the k·p expansion using combinatorial numbering:
%   - order=0: Constant term (index 1)
%   - order>0: Range of indices for that order's terms
%
% Algorithm:
%   Uses nchoosek to compute combinatorial ranges
%   Follows standard polynomial term ordering convention
%
% Example:
%   ord2_terms = orderlist(2); % Gets indices for all 2nd-order terms
%
% See also:
%   nchoosek
if order == 0
    Orderlist = 1;
else
    Orderlist = nchoosek(4+order-2,order-1)+1:nchoosek(4+order-1,order);
end
end
