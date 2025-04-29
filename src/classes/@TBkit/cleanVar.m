function SymListRow = cleanVar(SymListRow,Accuracy)
%CLEANVAR Simplify symbolic expressions by removing small terms
%
%   Syntax:
%       SymListRow = cleanVar(SymListRow,Accuracy)
%
%   Description:
%       Removes terms smaller than specified accuracy from symbolic
%       expressions and simplifies remaining terms.
%
%   Inputs:
%       SymListRow - Symbolic expression or array
%       Accuracy   - Threshold for term removal (default=1e-12)
%
%   Output:
%       SymListRow - Simplified symbolic expression
Accur  = round(-log(Accuracy)/log(10));
nSymListRow = length(SymListRow);
notChoose = true(ones(1,nSymListRow));
for i = 1:nSymListRow
    [coeffs_list,symvar_L] = coeffs(SymListRow(i));
    jChoose = false;
    sumCoeffsTmp = 0;
    for j = 1:length(coeffs_list)
        CoeffsTmp = coeffs(coeffs_list(j));
        try
            dCoeffsTmp = double(vpa(CoeffsTmp,Accur));
            if dCoeffsTmp >1e30
                warning('addtional numerical error!');
            end
        catch
            jChoose = true;
            break;
        end
        if ~isempty(dCoeffsTmp)
            sumCoeffsTmp = sumCoeffsTmp +abs(dCoeffsTmp);
        end
    end

    if ~jChoose
        if sumCoeffsTmp>Accuracy
            dcoeffs_list = roundn(double(coeffs_list),-Accur);
            jChoose = true;
            SymListRow(i) = sum(dcoeffs_list.*symvar_L);
            %disp(sumCoeffsTmp);
            %disp(vpa(SymListRow(i),6));
        end
    end
    notChoose(i) = ~jChoose;
end
SymListRow(notChoose) = sym(0);
end