function H_hr = charalize(H_hr)
for i = 1:H_hr.WAN_NUM
for j = 1:H_hr.WAN_NUM
if H_hr.elementL(i) == H_hr.elementL(j) && i~=j
H_hr.HnumL(i,j,:) = 0;
end
end
end
end
