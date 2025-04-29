function [H_hk,Sublist,Unique_term] =  unique(H_hk,seed,checklist,options)
%UNIQUE Simplify Hamiltonian coefficients by identifying unique terms
%
%   Syntax:
%       [H_hk,Sublist,Unique_term] = unique(H_hk,seed,checklist,options)
%
%   Description:
%       Identifies and replaces repeated coefficients in the Hamiltonian
%       with symbolic variables to simplify the expression.
%
%   Inputs:
%       H_hk      - HK object containing Hamiltonian
%       seed      - Base name for symbolic variables (default='eta')
%       checklist - List of coefficients to check for relations
%       options   - Struct with Accuracy and simplify fields
%
%   Outputs:
%       H_hk        - Modified HK object
%       Sublist     - Substitution list
%       Unique_term - List of unique terms found
arguments
    H_hk HK;
    seed char = 'eta';
    checklist =sym([-1,sqrt(3),-sqrt(3),2,-2,3,-3,6,-6]);
    options.Accuracy = 1e-6;
    options.simplify logical=true;
end
HcoeL_tmp_1 = H_hk.HcoeL(:);
HcoeL_tmp = [real(HcoeL_tmp_1);imag(HcoeL_tmp_1)];
if  options.simplify
    HcoeL_tmp = TBkit.cleanVar(HcoeL_tmp,options.Accuracy);
end
[Unique_term,~,ic] = unique(HcoeL_tmp);
SymVar = sym(seed,[length(Unique_term),1],'real');
SymVar(find(Unique_term==sym(0))) = sym(0);
for coeff_for_check = checklist
    [Lia,Locb] = ismember(expand(Unique_term),expand(coeff_for_check*Unique_term));
    Equationlist = HR.isolateAll((SymVar(Lia) == coeff_for_check*SymVar(Locb(Locb>0)) ));
    SymVar = subs(SymVar,lhs(Equationlist),rhs(Equationlist));
end
Sublist = SymVar == Unique_term;
HcoeL_tmp_1 = SymVar(ic(1:end/2)) ...
    +1i*SymVar(ic(end/2+1:end));
H_hk.HcoeL = reshape(HcoeL_tmp_1,H_hk.Basis_num,H_hk.Basis_num,H_hk.Kinds);
end
