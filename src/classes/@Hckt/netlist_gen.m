function [HcktObj,ScktnameL,AllPortStrL,ijkL] = netlist_gen(HcktObj,fid,options)
arguments
HcktObj Hckt;
fid = 1;
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
for t = 1:HcktObj.dim
dimMesh{t}  = (1:HcktObj.mesh(t)).';
end
if options.BinBin
PPdist = HcktObj.Port2PortDist_BinBin;
else
PPdist = HcktObj.Port2PortDist;
end
UsefulMesh = HcktObj.mesh(HcktObj.mesh>1);
Ntotmesh = prod(UsefulMesh);
ijkL = zeros(Ntotmesh,HcktObj.dim);
count = 0;
for t = find(HcktObj.mesh>1)
count = count +1;
if count >1
oldMesh = prod(UsefulMesh(1:count-1));
else
oldMesh = 1;
end
ijkL(:,t) = repmat(kron(dimMesh{t},ones(oldMesh,1)),...
[Ntotmesh/(HcktObj.mesh(t)*oldMesh),1]);
end
ScktnameL = repmat('X',[Ntotmesh,1]);
underscore_L=repmat('_',[Ntotmesh,1]);
nodepreL = repmat('n',[Ntotmesh,1]);
for i = 1:size(ijkL,2)
tmpcharL = park.num2strwithzero(ijkL(:,i));
ScktnameL = [ScktnameL,underscore_L,tmpcharL];
end
ScktnameL = string(ScktnameL);
portL{HcktObj.Nports} = ijkL;
AllPortStrL = [];
if strcmp(options.convention,'II')
for i = 1:HcktObj.Nports
ConnectionPort = PPdist(i);
portL{i} = portL{HcktObj.Nports} + HcktObj.Port2VectorDist(i);
portL{i}(portL{i}>HcktObj.mesh) = 1;
for j = 1:HcktObj.dim
portL{i}((portL{i}(:,j) ==0),j)= HcktObj.mesh(j);
end
tmpPortStrL= nodepreL;
for j = 1:HcktObj.dim
tmpcharL = park.num2strwithzero(portL{i}(:,j));
tmpPortStrL =[tmpPortStrL,underscore_L,tmpcharL];
end
tmpPortStrL = [tmpPortStrL,underscore_L,park.num2strwithzero(ConnectionPort*ones(Ntotmesh,1))];
PortStrL{i}  = tmpPortStrL;
AllPortStrL = [AllPortStrL,string(tmpPortStrL)];
end
else
for i = 1:HcktObj.Nports
ConnectionPort = i;
portL{i} = ijkL;
tmpPortStrL= nodepreL;
for j = 1:HcktObj.dim
tmpcharL = park.num2strwithzero(portL{i}(:,j));
tmpPortStrL =[tmpPortStrL,underscore_L,tmpcharL];
end
tmpPortStrL = [tmpPortStrL,underscore_L,park.num2strwithzero(ConnectionPort*ones(Ntotmesh,1))];
PortStrL{i}  = tmpPortStrL;
AllPortStrL = [AllPortStrL,string(tmpPortStrL)];
end
end
AllPortStrL = unique(AllPortStrL);
if strcmp(options.convention,'II')
basePriVname = 'PriV';
for i = 1:HcktObj.dim
basePriVname = [basePriVname,'_1'];
end
else
basePriVname = 'Pri';
end
PriVnameL = repmat(basePriVname,[Ntotmesh,1]);
for i = 1:Ntotmesh
fprintf(fid,"%s ",ScktnameL(i,:));
for n = 1:HcktObj.Nports
fprintf(fid,"%s ",PortStrL{n}(i,:));
end
fprintf(fid,"%s %s %s\n",'GND',PriVnameL(i,:),HcktObj.AddtionInformation);
end
end
