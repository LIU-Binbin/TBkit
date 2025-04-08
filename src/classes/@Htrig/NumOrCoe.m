function [num_label,coe_label,H_htrig] = NumOrCoe(H_htrig)
if isempty(H_htrig.num) ||  isempty(H_htrig.coe)
[num_label,coe_label] = NumOrCoe@TBkit(H_htrig);
H_htrig.num = num_label;
H_htrig.coe = coe_label;
else
num_label = H_htrig.num;
coe_label = H_htrig.coe;
end
end
