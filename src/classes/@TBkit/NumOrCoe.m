function [num_label,coe_label,H_hk] = NumOrCoe(H_hk)
if isempty(H_hk.num) ||  isempty(H_hk.coe)
[num_label,coe_label] = NumOrCoe@TBkit(H_hk);
H_hk.num = num_label;
H_hk.coe = coe_label;
else
num_label = H_hk.num;
coe_label = H_hk.coe;
end
end
