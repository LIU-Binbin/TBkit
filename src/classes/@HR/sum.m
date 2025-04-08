function H_hr = sum(H_hr_list)
H_hr = H_hr_list(1);
for i = 2:length(H_hr_list)
H_hr = H_hr + H_hr_list(i);
end
end
