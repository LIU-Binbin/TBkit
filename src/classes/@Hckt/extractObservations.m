function [VectorList,ObservationsMat,Vstruct,SpectrumX,TimeL] = extractObservations(simulation_result,Observations,options)
arguments
simulation_result struct;
Observations = 'v';
options.analysis {mustBeMember(options.analysis,{'tran','ac'})} = 'tran';
options.mode = 'simple';
options.x = true;
options.fft = true;
options.patten = "_"|".";
options.Vstrcut = false;
options.SelectNode = [];
end
chooseL = ones(1,length(simulation_result));
for i = 1:numel(simulation_result)
if ~strcmp(simulation_result(i).var_name(1),Observations)
chooseL(i) = false;
elseif strcmp(simulation_result(i).var_name ,[Observations,'(0)'])
chooseL(i) = false;
elseif ~options.x
if contains(simulation_result(i).var_name ,'x')
chooseL(i) = false;
end
else
end
end
switch  options.analysis
case 'tran'
TimeL =  simulation_result(1).val;
DurationTime = max(TimeL);
TimeListN = round(TimeL*(1/DurationTime)*length(TimeL));
SpectrumX = (1:length(TimeL))/(DurationTime);
simulation_result_Observations = simulation_result(logical(chooseL));
VectorList = [];
for i = 1:size(simulation_result_Observations,1)
tmpChar = simulation_result_Observations(i).var_name;
tmpChar(1:2) = [];
tmpChar(end) = [];
tmlStrL = split(tmpChar,options.patten).';
tmlStrL(1) = [];
try
VectorList = [VectorList;double([string(tmlStrL)])];
catch
VectorList = [VectorList;zeros(size(VectorList(1,:)))];
end
end
if isempty(options.SelectNode)
else
SelectL = VectorList(:,end) == options.SelectNode;
VectorList = VectorList(SelectL,:);
simulation_result_Observations = simulation_result_Observations(SelectL);
end
if options.Vstrcut
ObservationsMat = [];
for i = 1:size(VectorList,1)
Vstruct(i).R = VectorList(i,1:end-1);
Vstruct(i).port = VectorList(i,end);
Vstruct(i).val = simulation_result_Observations(i).val;
Vstruct(i).fftval = (nufft(simulation_result_Observations(i).val.',TimeListN));
if options.fft
ObservationsMat = [ObservationsMat;Vstruct(i).fftval];
else
ObservationsMat = [ObservationsMat;Vstruct(i).val.' ];
end
end
else
Vstruct = [];
ObservationsMat = zeros(size(VectorList,1),length(TimeL));
if options.fft
for i = 1:size(VectorList,1)
ObservationsMat(i,:) = (nufft(simulation_result_Observations(i).val.',TimeListN));
end
else
for i = 1:size(VectorList,1)
ObservationsMat(i,:) = simulation_result_Observations(i).val.' ;
end
end
end
case 'ac'
SpectrumX =  simulation_result(1).val;
TimeL = [];
simulation_result_Observations = simulation_result(logical(chooseL));
VectorList = [];
NPROBE = size(simulation_result_Observations,1);
IsRealL = ones(NPROBE,1);
for i = 1:NPROBE
tmpChar = simulation_result_Observations(i).var_name;
if tmpChar(2) == 'i'
IsRealL(i) = 1i;
end
tmpChar([1:3]) = [];
tmpChar(end) = [];
tmlStrL = split(tmpChar,options.patten).';
tmlStrL(1) = [];
try
VectorList = [VectorList;double([string(tmlStrL)])];
catch
VectorList = [VectorList;zeros(size(VectorList(1,:)))];
end
end
if isempty(options.SelectNode)
else
SelectL = VectorList(:,end) == options.SelectNode;
VectorList = VectorList(SelectL,:);
IsRealL = IsRealL(SelectL);
simulation_result_Observations = simulation_result_Observations(SelectL);
end
Vstruct = [];
ObservationsMat = zeros(size(VectorList,1),length(SpectrumX));
for i = 1:size(VectorList,1)
ObservationsMat(i,:) = IsRealL(i) * simulation_result_Observations(i).val.' ;
end
[VectorList,ObservationsMat] = HollowKnight.generalcontractrow2(VectorList,ObservationsMat);
otherwise
end
end
