function HcktObj = vectorize(HcktObj,fin_dir)
arguments
HcktObj
fin_dir = zeros(HcktObj.dim,1);
end
basePriVname = 'PriV';
for i = 1:HcktObj.dim
basePriVname = [basePriVname,'_1'];
end
count = 0;
checkvectorlist = HcktObj.dim2vectorL(HcktObj.dim);
checkvectorlist(1,:) = [];
count = count +1;
netlistTmp{count} = strcat(HcktObj.ScktL(1).type,num2str(count));
for i = 1:length(HcktObj.PortInCell{1})
netlistTmp{count} = [...
netlistTmp{count},...
string(HcktObj.PortInCell{1}(i)),...
];
end
for i = 1:length(HcktObj.PortInCell{1})
netlistTmp{count} = [...
netlistTmp{count},...
string((HcktObj.PortOutCell{1}(i))),...
];
end
if HcktObj.ScktL(1).ToGND
netlistTmp{count} = [netlistTmp{count},'TOGND'];
end
netlistTmp{count} = [netlistTmp{count},string(HcktObj.ScktL(1).name)];
for i = 2:size(HcktObj.ScktL,1)
tmpVector = HcktObj.Port2VectorDist(HcktObj.PortInCell{i});
[~,seq] = ismember(tmpVector,checkvectorlist,'rows');
if seq ~= 0
count = count +1;
netlistTmp{count} = strcat(HcktObj.ScktL(i).type,num2str(count));
for j = 1:length(HcktObj.PortInCell{i})
netlistTmp{count} = [...
netlistTmp{count},...
"n"+string(HcktObj.Port2PortDistForVectorise(HcktObj.PortOutCell{i}(j))),...
string(HcktObj.PortInCell{i}(j))];
end
if HcktObj.ScktL(i).ToGND
netlistTmp{count} = [netlistTmp{count},'TOGND'];
end
netlistTmp{count} = [netlistTmp{count},string(HcktObj.ScktL(i).name)];
end
end
HcktObj.VectorizeScktL = Subckt('X',basePriVname,...
["n"+([string(HcktObj.PortInCell{1})]),([string(HcktObj.PortOutCell{1})]),...
'TOGND'],...
HcktObj.HomeCell.Description,netlistTmp);
if sum(fin_dir>0)
AddingGNDScktL =sum(fin_dir>0);
end
end
