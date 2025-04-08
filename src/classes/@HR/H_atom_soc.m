function H_atom_soc = H_atom_soc(H_hr)
if isempty(H_hr.quantumL)
error('you should provide quantum number list');
end
if size(H_hr.quantumL,1) ~= H_hr.WAN_NUM
error('size quantum number list qrong');
end
if mod(size(H_hr.quantumL,1),2) ~= 0
error('size quantum number list odd, cant be spinful');
end
if sum(H_hr.quantumL(:,4)) ~= 0
error('spin up ~= spin dn');
end
H_atom_soc = sym(zeros(H_hr.WAN_NUM));
for i = 1:H_hr.WAN_NUM
l1 = H_hr.quantumL(i,2);
m1 = H_hr.quantumL(i,3);
s1 = H_hr.quantumL(i,4);
element1 = H_hr.elementL(i);
orb1 = H_hr.orbL(i,:);
for j =  1:H_hr.WAN_NUM
l2 = H_hr.quantumL(j,2);
m2 = H_hr.quantumL(j,3);
s2 = H_hr.quantumL(j,4);
element2 = H_hr.elementL(j);
orb2 = H_hr.orbL(j,:);
if isequal(orb1,orb2) && element1 == element2 && l1 == l2
H_atom_soc(i,j) = soc_term_gen(l1,l2,m1,m2,s1,s2,element1);
end
end
end
end
