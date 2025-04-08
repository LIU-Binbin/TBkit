function [H_hr,Sublist,Unique_term] = unique(H_hr,seed,checklist,Accuracy)
arguments
H_hr HR;
seed char = 'gamma';
checklist =sym([sqrt(3)]);
Accuracy = (1e-6);
end
if strcmp(H_hr.Type,'list')
HcoeL_tmp = [real(H_hr.HcoeL);imag(H_hr.HcoeL)];
[CoeffsList,factor_list_2] = TBkit.factorAll(HcoeL_tmp);
[Unique_term,~,ic] = unique(factor_list_2);
SymVar = sym(seed,[length(Unique_term),1],'real');
SymVar(find(Unique_term==sym(0))) = sym(0);
for coeff_for_check = checklist
[Lia,Locb] = ismember(expand(Unique_term),expand(coeff_for_check*Unique_term));
Equationlist = HR.isolateAll((SymVar(Lia) == coeff_for_check*SymVar(Locb(Locb>0)) ));
SymVar = subs(SymVar,lhs(Equationlist),rhs(Equationlist));
end
Sublist = SymVar == Unique_term;
H_hr.HcoeL = SymVar(ic(1:end/2)) .* CoeffsList(1:end/2) ...
+1i*SymVar(ic(end/2+1:end)) .* CoeffsList(end/2+1:end) ;
end
end
