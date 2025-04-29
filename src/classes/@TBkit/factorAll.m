function [factor_list_1,factor_list_2] = factorAll(SymListRow)
%FACTORALL Factorize symbolic expressions into two components
%   [f1, f2] = factorAll(SymListRow) decomposes each symbolic expression
%   into two factors when possible.
%
%   Input:
%       SymListRow - Row vector of symbolic expressions
%
%   Outputs:
%       factor_list_1 - First factors (defaults to 1)
%       factor_list_2 - Second factors (contains remaining terms)
%
%   Example:
%       [f1, f2] = factorAll([sym('x*y'), sym('2*z')]);
    factor_list_1 = sym(ones(size(SymListRow)));
    factor_list_2 = sym(zeros(size(SymListRow)));
    
    for i = 1:length(SymListRow)
        F = factor(SymListRow(i));
        if length(F)==2
            factor_list_1(i) = F(1);
            factor_list_2(i)  = F(2);
        else
            factor_list_2(i)  = F(1);
        end
    end
end