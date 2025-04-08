function H_hk = hermitize(H_hk)
if isequal(zeros(size(H_hk.HnumL)),H_hk.HnumL)
num_label = false;
else
num_label = true;
end
if isequal(sym(zeros(size(H_hk.HcoeL))),H_hk.HcoeL)
coe_label = false;
else
coe_label = true;
end
H_hk_bk = H_hk';
if coe_label
H_hk = (H_hk + H_hk_bk)/2;
end
if num_label
H_hk.HnumL = (H_hk_bk.HnumL + H_hk.HnumL )/2;
end
end
