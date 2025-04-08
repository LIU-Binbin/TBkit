function H_htrig = find_HsymL_trig_bk(H_htrig,coeff_trig)
VarUsing = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
for j =1:length(coeff_trig)
coeff_trig_str = string(sym('A')*coeff_trig(j));
if isequal(H_htrig.seeds,string(VarUsing)) && strcmp(H_htrig.Type,'sincos')
pat1 = "*"+("sin"|"cos");
pat2 = "^"+digitsPattern(1) + "*"+("sin"|"cos");
Type = 'sincos';
elseif sum(contains(H_htrig.seeds,"Sigma_x"|"Sigma_y"|"Sigma_z"|"Sigma_w")) && strcmp(H_htrig.Type,'slab')
pat1 = "*"+"Sigma_";
Type = 'slab';
elseif isequal(H_htrig.seeds,string(VarUsing)) && strcmp(H_htrig.Type,'exp')
pat1 = "*"+"exp";
Type = 'exp';
else
end
if strcmp(Type,'sincos')
coeff_trig_str_list = split(coeff_trig_str,[pat1,pat2]);
for i = 2:numel(coeff_trig_str_list)
if coeff_trig_str_list(i) ~= ""
sym1 = str2sym(strcat('cos',coeff_trig_str_list(i)));
if ~ismember(sym1,H_htrig.HsymL_trig_bk)
sym2 = str2sym(strcat('sin',coeff_trig_str_list(i)));
H_htrig.HsymL_trig_bk = [H_htrig.HsymL_trig_bk,sym1,sym2];
end
end
end
elseif strcmp(Type,'slab')
coeff_trig_str_list = split(coeff_trig_str,[pat1]);
for i = 2:numel(coeff_trig_str_list)
if coeff_trig_str_list(i) ~= ""
if ~contains(string(H_htrig.HsymL_trig_bk),coeff_trig_str_list(i))
sym1 = str2sym(strcat('Sigma_',coeff_trig_str_list(i)));
H_htrig.HsymL_trig_bk = [H_htrig.HsymL_trig_bk,sym1];
end
end
end
elseif strcmp(Type,'exp')
coeff_trig_str_list = split(coeff_trig_str,[pat1]);
for i = 2:numel(coeff_trig_str_list)
if coeff_trig_str_list(i) ~= ""
if ~contains(string(H_htrig.HsymL_trig_bk),coeff_trig_str_list(i))
sym1 = str2sym(strcat('exp',coeff_trig_str_list(i)));
sym2 = str2sym(strcat('exp(-',coeff_trig_str_list(i),")"));
H_htrig.HsymL_trig_bk = [H_htrig.HsymL_trig_bk,sym1,sym2];
end
end
end
else
end
end
end
