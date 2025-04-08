function [coeff_trig,symvar_list_trig,H_htrig] = split_sym_eq(H_htrig,symbolic_polynomial)
switch H_htrig.Type
case {'exp'}
k_symbol = combine(rewrite(symbolic_polynomial,'exp'));
case {'mat','list'}
k_symbol = expand(simplify(rewrite(symbolic_polynomial,'exp'),'IgnoreAnalyticConstraints',true));
case 'sincos'
k_symbol = expand(simplify(rewrite(symbolic_polynomial,'sincos'),'IgnoreAnalyticConstraints',true));
otherwise
k_symbol = expand(rewrite(simplify(symbolic_polynomial),'sincos'));
end
k_symbol_str = string(k_symbol);
symvar_list = symvar(k_symbol);
for i = 1:length(symvar_list)
str_tmp = string(symvar_list(i));
switch str_tmp
case {'k_x','k_X','K_X','K_x','kx','kX','KX','Kx'}
k_symbol_str = strrep(k_symbol_str,str_tmp,'k_x');
case {'k_y','k_Y','K_y','K_Y','ky','kY','Ky','KY'}
k_symbol_str = strrep(k_symbol_str,str_tmp,'k_y');
case {'k_z','k_Z','K_z','K_Z','kz','kZ','Kz','KZ'}
k_symbol_str = strrep(k_symbol_str,str_tmp,'k_z');
case {'k_w','k_W','K_w','K_W','kw','kW','Kw','KW'}
k_symbol_str = strrep(k_symbol_str,str_tmp,'k_w');
end
end
switch H_htrig.Type
case {'exp','mat','list'}
k_symbol = combine(str2sym(k_symbol_str),'exp');
otherwise
k_symbol = expand(combine(str2sym(k_symbol_str),'sincos'));
end
k_symbol_children1 = children(k_symbol);
if isequal(simplify(fold(@mtimes,[k_symbol_children1{:}])),k_symbol)
k_symbol_children{1} = k_symbol;
elseif isequal(simplify(fold(@plus,[k_symbol_children1{:}])),k_symbol)
k_symbol_children = k_symbol_children1;
else
k_symbol_children{1} = k_symbol;
end
switch H_htrig.Type
case {'mat','list'}
nSon = length(k_symbol_children);
symvar_list_trig = zeros(nSon,3,'sym');
coeff_trig = zeros(nSon,1,'sym');
if k_symbol == sym(0)
else
for k = 1:nSon
[TmpCoe,TmpVar] = coeffs(k_symbol_children{k});
ActualVar = simplify(log(TmpVar),'IgnoreAnalyticConstraints',true)/1i;
[ActualCoe,ktype] = coeffs(ActualVar,H_htrig.seedsvar);
if ismember(sym(1),ktype)
n1 = find(sym(1)==ktype);
TmpCoe = simplify(exp(ActualCoe(n1)*(1i)),'IgnoreAnalyticConstraints',true)*TmpCoe;
ActualCoe(n1) = [];
ktype(n1) = [];
end
HsymL_coeLtmp = zeros(1,3,'sym');
for i = 1:3
ik = find(H_htrig.seedsvar(i) == ktype);
if ~isempty(ik)
HsymL_coeLtmp(i) = ActualCoe(ik);
end
end
symvar_list_trig(k,:) = HsymL_coeLtmp;
coeff_trig(k) = TmpCoe;
end
end
if strcmp(H_htrig.Type,'mat')
H_htrig = H_htrig.add_empty_one(symvar_list_trig);
end
return;
otherwise
for k = 1:length(k_symbol_children)
try
[coeff_trig,~] = coeffs(k_symbol_children{k},H_htrig.HsymL_trig_bk);
catch
coeff_trig = k_symbol_children{k};
end
for i = 1:numel(coeff_trig)
tmp_label = contains(string(coeff_trig),H_htrig.seeds);
if sum(tmp_label)
[~,coeff_trig_list] = coeffs(coeff_trig(i));
for j = 1:numel(coeff_trig_list)
tmp_label2 = contains(string(coeff_trig_list(j)),H_htrig.seeds);
if sum(tmp_label2)
H_htrig = H_htrig.find_HsymL_trig_bk(coeff_trig_list(j));
end
end
end
end
end
end
coeff_trig = sym([]);
symvar_list_trig= sym([]);
for k = 1:length(k_symbol_children)
[coeff_trig_tmp,symvar_list_trig_tmp] = coeffs(k_symbol_children{k},H_htrig.HsymL_trig_bk);
coeff_trig = [coeff_trig,coeff_trig_tmp];
symvar_list_trig = [symvar_list_trig,symvar_list_trig_tmp];
end
end
