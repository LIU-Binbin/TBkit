function HcktObj = autohermi(HcktObj)
arguments
HcktObj Hckt;
end
HcktObj_tmp = HcktObj;
for i = 2:HcktObj.NRPTS
vector_tmp = HcktObj.vectorAll(HcktObj.vectorL(i),:);
portin_tmp = HcktObj.PortInCell{i,1};
portout_tmp = HcktObj.PortOutCell{i,1};
vector_tmp_oppo = -vector_tmp;
portin_tmp_oppo = portout_tmp;
portout_tmp_oppo = portin_tmp;
[~,vector_tmp_oppo_label]= ismember(vector_tmp_oppo,HcktObj.vectorAll,'rows');
if vector_tmp_oppo_label == 0
fprintf('The opposite vector hamilton does not exist, build it : %s!\n', ...
string(HcktObj.ScktL(i).name));
HcktObj_tmp = HcktObj_tmp.set_hop(vector_tmp_oppo,HcktObj.ScktL(i),portin_tmp_oppo,portout_tmp_oppo);
continue;
end
vector1 = vector_tmp_oppo_label == HcktObj.vectorL;
i1 = park.checkcell(portin_tmp_oppo,HcktObj.PortInCell);
j1 = park.checkcell(portout_tmp_oppo,HcktObj.PortInCell);
j = find(i1&j1&vector1);
if i == j
end
if j == 0
fprintf('The opposite vector hamilton does not exist, build it : %s!\n', ...
string(HcktObj.ScktL(i).name));
HcktObj_tmp = HcktObj_tmp.set_hop(vector_tmp_oppo,HcktObj.HcktObj.ScktL(i),portin_tmp_oppo,portout_tmp_oppo);
continue;
end
if ~(HcktObj.ScktL(i) == HcktObj.ScktL(j))
disp([vector_tmp vector_tmp_oppo]);
fprintf('The opposite vector Subckt does not exist, build it with : %s vs %s!\n',HcktObj.ScktL(i).name,HcktObj.ScktL(j).name);
end
end
HcktObj = HcktObj_tmp;
end
