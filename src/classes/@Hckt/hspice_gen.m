function [HcktObj,Basename] = hspice_gen(HcktObj,filename,ops,options,options2,options3,options4,options5)
arguments
HcktObj Hckt;
filename = "";
ops.checking = true;
ops.commute = true;
options.mesh = repmat(10,[HcktObj.dim,1]);
options.fin_dir = zeros(HcktObj.dim,1);
options2.probenode {mustBeMember(options2.probenode,{'Allnode','Selectnode'})}= 'Allnode';
options2.prefix = "+ ";
options2.father = "";
options2.nodelist = [];
options3.fft = false;
options3.libmode = false;
options4.ExciteNodes = 1;
options4.mode {mustBeMember(options4.mode,{'general','vectorized','BinBin'})}= 'general';
options5.analysis {mustBeMember(options5.analysis,{'tran','ac'})} = 'tran';
end
optionscell = namedargs2cell(options);
options2cell = namedargs2cell(options2);
options3cell = namedargs2cell(options3);
options5cell = namedargs2cell(options5);
NodeStr  = 'n';
Basename = HcktObj.title;
for  i = 1:length(options.mesh)
NodeStr = [NodeStr,'_',num2str(options.mesh(i))];
Basename = [Basename,'_',num2str(options.mesh(i))];
end
if size(options4.ExciteNodes,2) == 1
for i = 1:size(options4.ExciteNodes,1)
NodeStrList{i} = [NodeStr,'_',num2str(options4.ExciteNodes(i))];
end
elseif size(options4.ExciteNodes,2) >1
for i = 1:size(options4.ExciteNodes,1)
NodeStr  = 'n';
for j = 1: size(options4.ExciteNodes,2)-1
NodeStr = [NodeStr,'_',num2str(options4.ExciteNodes(i,j))];
end
NodeStr = [NodeStr,'_',num2str(options4.ExciteNodes(i,j+1))];
NodeStrList{i} = NodeStr;
end
end
ICSTRING = ['.ic v(',NodeStrList{1},') = 1'];
switch HcktObj.magnitude
case 'p'
IpulseSTRING = ['Ipulse ',NodeStrList{1},' GND PU 0 1 5n 5n 50u'];
case 'u'
IpulseSTRING = ['Ipulse ',NodeStrList{1},' GND PU 0 1 5n 5n 50m'];
case 'm'
IpulseSTRING = ['Ipulse ',NodeStrList{1},' GND PU 0 1 5n 5n 500m'];
case 'n'
IpulseSTRING = ['Ipulse ',NodeStrList{1},' GND PU 0 1 5n 5n 500u'];
end
VacSTRING = "Vac source GND AC 1 0";
VacSTRING = [VacSTRING;string(['R_for_ac ',NodeStrList{1},' source 100'])];
for i = 2:length(NodeStrList)
ICSTRING = [ICSTRING;string(['.ic v(',NodeStrList{i},') = 1'])];
switch HcktObj.magnitude
case 'p'
IpulseSTRING = [IpulseSTRING;string(['Ipulse',num2str(i),' ',NodeStrList{i},' GND PU 0 1 5n 5n 50u'])];
case 'u'
IpulseSTRING = [IpulseSTRING;string(['Ipulse',num2str(i),' ',NodeStrList{i},' GND PU 0 1 5n 5n 50m'])];
case 'm'
IpulseSTRING = [IpulseSTRING;string(['Ipulse',num2str(i),' ',NodeStrList{i},' GND PU 0 1 5n 5n 500m'])];
case 'n'
IpulseSTRING = [IpulseSTRING;string(['Ipulse',num2str(i),' ',NodeStrList{i},' GND PU 0 1 5n 5n 50u'])];
end
VacSTRING = [VacSTRING;"Vac"+num2str(i)+" source"+ num2str(i)+" GND AC 1 0"];
VacSTRING = [VacSTRING;string(['R_for_ac',num2str(i),' ',NodeStrList{i},' source',num2str(i),' 100'])];
end
if strcmp(HcktObj.Vac,"")
HcktObj.Vac = VacSTRING;
end
if strcmp(HcktObj.IC,"")
HcktObj.IC = string(ICSTRING);
end
if strcmp(HcktObj.Ipulse,"")
HcktObj.Ipulse = string(IpulseSTRING);
end
if strcmp(filename,"")
switch options5.analysis
case 'ac'
filename = strcat(Basename,'_AC.sp');
otherwise
filename = strcat(Basename,'.sp');
end
end
switch options4.mode
case 'vectorized'
HcktObj = hspice_gen_vectorized(HcktObj,filename,optionscell{:},options2cell{:},options3cell{:});
case 'BinBin'
HcktObj = hspice_gen_vectorized(HcktObj,filename,optionscell{:},options2cell{:},options3cell{:},'BinBin',true);
case 'general'
if ops.checking
HcktObj = HcktObj.autohermi();
end
if ops.commute
HcktObj = HcktObj.half(ops.commute);
else
HcktObj = HcktObj.half();
end
HcktObj = hspice_gen_general(HcktObj,filename,optionscell{:},options2cell{:},options3cell{:},options5cell{:});
end
end
