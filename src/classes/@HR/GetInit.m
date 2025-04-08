function [VarInit,EQL2] =GetInit(H_hr,H_hr2,vectorL)
H_hr = H_hr.rewrite();
H_hr = H_hr.simplify();
H_hr2 =H_hr2.rewrite();
if nargin < 3
[ia,ic] = ismember(H_hr.vectorL,H_hr2.vectorL,'row');
HcoeLtmp = H_hr.HcoeL(ia);
HnumLtmp = H_hr2.HnumL(ic);
else
if size(vectorL,1) == 1
ic = all((vectorL) == H_hr.vectorL(:,1:H_hr.Dim),2);
HcoeLtmp = H_hr.HcoeL(ic);
vector1 = H_hr.vectorL(ic,:);
ic = all((vectorL) == H_hr2.vectorL(:,1:H_hr.Dim),2);
HnumLtmp = H_hr2.HnumL(ic);
vector2 = H_hr2.vectorL(ic,:);
[ia,ic] = ismember(vector1,vector2,'row');
HcoeLtmp = HcoeLtmp(ia);
HnumLtmp = HnumLtmp(ic);
end
end
[Unique_term,ia,ic] = unique(HcoeLtmp);
for i  = 1:length(ia)
EQL(i,:) = Unique_term(i) == max(HnumLtmp(ic == i));
end
[Unique_term2,ia,ic] = unique(lhs(EQL));
HcoeLtmp2 = rhs(EQL);
for i  = 1:length(ia)
EQL2(i,:) = Unique_term2(i) == max(HcoeLtmp2(ic == i));
end
VarInit = solve(EQL2,H_hr.symvar_list);
EQL2 = vpa(EQL2);
end
