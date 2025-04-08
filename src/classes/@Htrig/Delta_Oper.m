function Delta_Oper = Delta_Oper(Oper_list)
tmpeq_all = "1";
syms i1 i2 i3 j1 j2 j3 integer;
tmpeqstr_list  =["(i1 == j1)","(i2 == j2)","(i3 == j3)"];
for i = 1:length(Oper_list)
tmpeq ="";
strtmp = strsplit(string(Oper_list(i)),{'_'},'CollapseDelimiters',true);
switch strtmp{2}
case 'x'
tmpeq = tmpeq+'i1 == (j1';
case 'y'
tmpeq = tmpeq+'i2 == (j2';
case 'z'
tmpeq = tmpeq+'i3 == (j3';
otherwise
end
if strcmp(strtmp{3},'N')
tmpeq = tmpeq + ")";
else
tmpeq = tmpeq + " + "+strtmp{3}+")";
end
tmpeq = "(" +tmpeq+")";
switch strtmp{2}
case 'x'
tmpeqstr_list(1) = tmpeq;
case 'y'
tmpeqstr_list(2) = tmpeq;
case 'z'
tmpeqstr_list(3) = tmpeq;
otherwise
end
end
tmpeq_all = tmpeq_all + "&" +...
tmpeqstr_list(1)+'&'+...
tmpeqstr_list(2)+'&'+...
tmpeqstr_list(3);
Delta_Oper = matlabFunction(str2sym(tmpeq_all),'Vars',[i1 i2 i3 j1 j2 j3]);
end
