function H_htrig2 = subs(H_htrig,varargin)
HsymL_trig_tmp = H_htrig.HsymL_trig;
switch length(varargin)
case 0
HsymL_trig_tmp = expand(simplify(subs(HsymL_trig_tmp)));
case 1
HsymL_trig_tmp = expand(simplify(subs(HsymL_trig_tmp,varargin{1})));
case 2
HsymL_trig_tmp = expand(simplify(subs(HsymL_trig_tmp,varargin{1},varargin{2})));
end
H_htrig2 = H_htrig;
H_htrig2.HcoeL = sym([]);
H_htrig2.HnumL = [];
H_htrig2.HsymL_trig = sym([]);
count = 0;
for i = 1:H_htrig.Kinds
[coeff_trig,symvar_list_trig,H_htrig2] = split_sym_eq(H_htrig2,HsymL_trig_tmp(i));
for j =1:numel(coeff_trig)
count = count+1;
k_cell{count} = symvar_list_trig(j);
mat_cell{count} = H_htrig.HcoeL(:,:,i);
Var_cell{count} = coeff_trig(j);
end
end
H_htrig2 = H_htrig2.setup(Var_cell,k_cell,mat_cell);
end
