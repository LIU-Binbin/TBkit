function Factorlist_parity = factorlist_parity(Degree)
Nterm = nchoosek(4+Degree-1,Degree);
Factorlist_parity = zeros(1,Nterm);
for i =0:Degree
if mod(i, 2) == 0
Factorlist_parity(HK.orderlist(i)) = 1;
else
Factorlist_parity(HK.orderlist(i)) = -1;
end
end
end
