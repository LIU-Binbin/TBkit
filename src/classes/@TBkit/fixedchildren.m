function [ChirdrenCell,Type] = fixedchildren(SymVar,mode)
%FIXEDCHILDREN Decompose symbolic expression with controlled behavior
%
%   Syntax:
%       [ChirdrenCell,Type] = fixedchildren(SymVar,mode)
%
%   Description:
%       Recursively decomposes symbolic expressions while preserving
%       mathematical structure (sums, products) according to specified mode.
%
%   Inputs:
%       SymVar - Symbolic expression to decompose
%       mode   - Decomposition mode ('exp_inner' or default)
%
%   Outputs:
%       ChirdrenCell - Cell array of decomposed components
%       Type         - Structure type ('sum', 'prod', or 'inner')
%
%   See also: combine, simplify
arguments
    SymVar sym;
    mode = 'exp_inner';
end
if strcmp(mode,'exp_inner')
    subexpr = children(combine(SymVar,'exp'));
    if isequal(simplify(fold(@plus,subexpr)-SymVar),sym(0))
        ChirdrenCell = subexpr;
        Type = 'sum';
        return;
    end
    if isequal(simplify(fold(@times,subexpr)-SymVar),sym(0))
        ChirdrenCell{1} = SymVar;
        Type = 'prod';
        return;
    end
    ChirdrenCell = subexpr;
    Type = 'inner';
else
    subexpr = children(SymVar);
    if isequal(simplify(fold(@plus,subexpr)-SymVar),sym(0))
        ChirdrenCell = subexpr;
        Type = 'sum';
        return;
    end
    if isequal(simplify(fold(@times,subexpr)-SymVar),sym(0))
        ChirdrenCell{1} = SymVar;
        Type = 'prod';
        return;
    end
    ChirdrenCell{1} = SymVar;
    Type = 'inner';
end
end