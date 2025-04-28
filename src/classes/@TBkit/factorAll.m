function [factor_list_1,factor_list_2] = factorAll(SymListRow)
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