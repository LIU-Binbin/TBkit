function H_hk = applyR(H_hk,R)
arguments
H_hk HK;
R ;
end
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
matcell = H_hk.matgen(R);
if num_label
for i = 0:H_hk.Degree
Orderlist = HK.orderlist(i);
H_hk.HnumL(:,:,Orderlist) = HK.matrixtimespage(matcell{i+1},H_hk.HnumL(:,:,Orderlist));
end
end
if coe_label
for i = 0:H_hk.Degree
Orderlist = HK.orderlist(i);
H_hk.HcoeL(:,:,Orderlist) = HK.matrixtimespage(matcell{i+1},H_hk.HcoeL(:,:,Orderlist));
end
end
end
