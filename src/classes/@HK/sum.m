function H_hk = sum(H_hk_list)
H_hk = H_hk_list(1);
for i = 2:length(H_hk_list)
H_hk = H_hk + H_hk_list(i);
end
end
