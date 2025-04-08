function [simulation_result]=read_hspice_tr(filename,options)
arguments
filename = "Hckt.tr0";
options.fast = true;
options.filenameInformation = '';
options.ErrorMode = 0;
end
if options.fast
[~,~,ext] = fileparts(filename);
file_extension=ext(1:3);
if strcmp(options.filenameInformation,'')
filenameInformation = char(filename);
filenameInformation = [filenameInformation,'.info'];
else
filenameInformation = options.filenameInformation;
end
fileID = fopen(filenameInformation);
content_file = textscan(fileID,'
fclose(fileID);
content_str=string(content_file{1,1});
data_start_ind=find(content_str==('$&
FirstUsefulStrLogicalIndex = find(isnan(str2double(content_str(1:data_start_ind-2))),1);
FirstUsefulLabel = (FirstUsefulStrLogicalIndex);
content_str(1:FirstUsefulLabel-1) = [];
num_var=FirstUsefulLabel-1;
data_start_ind = data_start_ind +1 -FirstUsefulLabel;
simulation_result = Hckt.save_signal_names(num_var,data_start_ind,content_str,file_extension);
if exist(filename+".data",'file')
fileID = fopen(filename+".data");
elseif exist(filename,'file')
fileID = fopen(filename);
end
formatSpec = '
dataArray = textscan(fileID, formatSpec);
checkEndRow = 6;
for i = 1:6
if isnan(dataArray{i}(end))
checkEndRow = i-1;
break;
end
end
iTOT = length(dataArray{1});
TOTnum = iTOT*6;
DATA = zeros(iTOT,6);
count = 1;
for i = 1:6
if i < checkEndRow
DATA(1:iTOT,i) = (dataArray{i}(1:end));
count = count+iTOT;
else
DATA(1:(iTOT-1),i) = (dataArray{i}(1:end-1));
count = count+iTOT-1;
end
end
DATA = DATA.';
DATA = DATA(:);
DATA(count:end)= [];
if  options.ErrorMode
DATA(1) = [];
end
NStruct= length(simulation_result);
DATA = reshape(DATA,NStruct,numel(DATA)/NStruct).';
for i = 1:length(simulation_result)
simulation_result(i).val = DATA(:,i);
end
fclose(fileID);
else
try
if  ~isunix
win_matlab.hspiceACTR(filename);
else
end
options.fast = true;
optionsCell = namedargs2cell(options);
[simulation_result]= Hckt.read_hspice_tr(filename,optionsCell{:});
catch
[simulation_result]= Hckt.read_hspice_tr_sw_ac(filename);
end
end
end
