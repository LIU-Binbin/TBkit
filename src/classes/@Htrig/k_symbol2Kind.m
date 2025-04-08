function Kind = k_symbol2Kind(H_htrig,k_symbol)
switch H_htrig.Type
case {'exp','sincos'}
str_2_compare = string(k_symbol);
Kind = find(string(H_htrig.HsymL_trig) == str_2_compare);
case {'mat','list'}
if isa(k_symbol,'sym')
[~,Kind]=ismember(k_symbol,H_htrig.HsymL_coeL,'rows');
else
[~,Kind]=ismember(k_symbol,H_htrig.HsymL_numL,'rows');
end
otherwise
str_2_compare = string(k_symbol);
Kind = find(string(H_htrig.HsymL_trig) == str_2_compare);
end
end
