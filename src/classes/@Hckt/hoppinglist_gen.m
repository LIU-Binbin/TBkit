function [ScktnameL] = hoppinglist_gen(HcktObj,fid,seq,ijkL,ScktnameL,options)
arguments
HcktObj Hckt;
fid = 1;
seq = 2;
ijkL = [];
ScktnameL = [];
options.BinBin = false;
options.mesh = repmat(10,[HcktObj.dim,1]);
options.fin_dir = zeros(HcktObj.dim,1);
options.convention = 'I';
end
switch class(fid)
case {'char','string'}
fid = fopen(fid,'w');
case 'double'
end
if isempty(HcktObj.mesh)
HcktObj.mesh = options.mesh ;
end
if length(HcktObj.mesh) ~= HcktObj.dim
error('wrong input: mesh! The dim of HcktObj is
end
if isempty(HcktObj.fin_dir)
HcktObj.fin_dir = options.fin_dir ;
end
if sum(HcktObj.fin_dir)
fin_dir_mode=true;
else
fin_dir_mode=false;
end
UsefulMesh = HcktObj.mesh(HcktObj.mesh>1);
Ntotmesh = prod(UsefulMesh);
ScktnameL = strrep(ScktnameL,'X',...
[HcktObj.ScktL(seq).type,num2str(seq)]);
underscore_L=repmat('_',[Ntotmesh,1]);
nodepreL = repmat('n',[Ntotmesh,1]);
NPortInCell = length(HcktObj.PortInCell{seq});
NPortOutCell = length(HcktObj.PortOutCell{seq});
NPort = NPortInCell+NPortOutCell;
portL{NPort} = ijkL;
GeneralRvector =  HcktObj.vectorAll(HcktObj.vectorL(seq),:);
if strcmp(options.convention,'I')
for i = 1:NPortInCell
ConnectionPort{i} = HcktObj.PortInCell{seq}(i);
portL{i} = ijkL;
end
if fin_dir_mode
DupliteportL = portL;
DupliteportLabel = false(length(ScktnameL),1);
exceed_logical{NPort} = false(length(ScktnameL),1);
end
for i = (NPortInCell+1):NPort
ConnectionPort{i} = HcktObj.PortOutCell{seq}(i-NPortInCell);
portL{i} = ijkL + GeneralRvector;
if fin_dir_mode
DupliteportL{i} = portL{i};
end
exceed_logical{i} = false(length(ScktnameL),1);
for j = 1:HcktObj.dim
exceed_r_logical = portL{i}(:,j)>HcktObj.mesh(j);
exceed_r_label = find(exceed_r_logical);
portL{i}(exceed_r_label,j) = mod(portL{i}(exceed_r_label,j),HcktObj.mesh(j));
if fin_dir_mode
if HcktObj.fin_dir(j)
exceed_logical{i} = exceed_logical{i} + exceed_r_logical;
DupliteportL{i}(exceed_r_label,j) = mod(portL{i}(exceed_r_label,j),HcktObj.mesh(j));
DupliteportL{i-NPortInCell}(exceed_r_label,j) = 0;
else
DupliteportL{i}(exceed_r_label,j) = portL{i}(exceed_r_label,j);
end
end
end
for j = 1:HcktObj.dim
exceed_l_logical = (portL{i}(:,j)-1) < 0 ;
exceed_l_label = find(exceed_l_logical);
portL{i}(exceed_l_label,j)= mod(portL{i}(exceed_l_label,j)-1,HcktObj.mesh(j));
if fin_dir_mode
if HcktObj.fin_dir(j)
exceed_logical{i} = exceed_logical{i} + exceed_l_logical;
DupliteportL{i}(exceed_l_label,j) = mod(portL{i}(exceed_l_label,j)-1,HcktObj.mesh(j));
DupliteportL{i-NPortInCell}(exceed_l_label,j) = 0;
else
DupliteportL{i}(exceed_l_label,j) = portL{i}(exceed_l_label,j);
end
end
end
if fin_dir_mode
portL{i}(logical(exceed_logical{i}) ,j) = 0 ;
DupliteportLabel = DupliteportLabel + exceed_logical{i};
end
end
for i = 1:NPort
tmpPortStrL= nodepreL;
for j = 1:HcktObj.dim
tmpcharL = park.num2strwithzero(portL{i}(:,j));
tmpPortStrL =[tmpPortStrL,underscore_L,tmpcharL];
end
tmpPortStrL = [tmpPortStrL,underscore_L,park.num2strwithzero(ConnectionPort{i}*ones(Ntotmesh,1))];
PortStrL{i}  = string(tmpPortStrL);
if i > NPortInCell
PortStrL{i}(logical(exceed_logical{i})) = "GND";
else
end
end
if fin_dir_mode
DupliteportLabel = logical(DupliteportLabel);
AddtionalScktnameL = ScktnameL(DupliteportLabel,:);
Naddtional = size(AddtionalScktnameL,1);
AddtionalScktnameL = strcat(AddtionalScktnameL,repmat('_extra',[Naddtional,1]));
AddtionalnodepreL = repmat('n',[Naddtional,1]);
Addtionalunderscore_L =repmat('_',[Naddtional,1]);
for i = 1:NPort
tmpPortStrL= AddtionalnodepreL;
for j = 1:HcktObj.dim
tmpcharL = park.num2strwithzero([DupliteportL{i}(DupliteportLabel,j);portL{i}(:,j)]);
tmpcharL = tmpcharL(1:Naddtional,:);
tmpPortStrL =[tmpPortStrL,Addtionalunderscore_L,tmpcharL];
end
ConnectionPortStrL = park.num2strwithzero(ConnectionPort{i}*ones(Naddtional,1));
tmpPortStrL = [tmpPortStrL,Addtionalunderscore_L,ConnectionPortStrL];
AddtionalPortStrL{i}  = string(tmpPortStrL);
if i <= NPortInCell
AddtionalPortStrL{i}(1:Naddtional) = "GND";
else
end
end
end
else
end
if strcmp(HcktObj.ScktL(seq).type,'X')
basePriVname = strcat("GND ",HcktObj.ScktL(seq).name);
else
basePriVname = strcat(HcktObj.ScktL(seq).name,HcktObj.ScktL(seq).description);
end
PriVnameL = repmat(basePriVname,[Ntotmesh,1]);
if fin_dir_mode
AddtionalPriVnameL = repmat(basePriVname,[Naddtional,1]);
end
for i = 1:Ntotmesh
fprintf(fid,"%s ",ScktnameL(i,:));
for n = 1:NPort
fprintf(fid,"%s ",PortStrL{n}(i,:));
end
fprintf(fid,"%s %s\n",PriVnameL(i,:),HcktObj.DescriptionL{seq});
end
if sum(HcktObj.fin_dir)
for i = 1:Naddtional
fprintf(fid,"%s ",AddtionalScktnameL(i,:));
for n = 1:NPort
fprintf(fid,"%s ",AddtionalPortStrL{n}(i,:));
end
fprintf(fid,"%s %s\n",AddtionalPriVnameL(i,:),HcktObj.DescriptionL{seq});
end
end
end
