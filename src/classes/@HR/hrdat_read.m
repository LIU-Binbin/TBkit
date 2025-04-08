function [dataArray,NRPT_list,NRPTS,NUM_WAN]=hrdat_read(filename)
if nargin < 1
filename = 'wannier90_hr.dat';
end
delimiter = ' ';
startRow = 2;
endRow = 3;
formatSpec = '
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter,...
'MultipleDelimsAsOne', true, 'TextType', 'string',...
'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
NBR = dataArray{:, 1};
NUM_WAN = NBR(1);
NRPTS = NBR(2);
NRPTS_num1=fix(double(NRPTS)/15);
NRPTS_num2=mod(NRPTS,15);
fclose(fileID);
fileID = fopen(filename,'r');
startRow = 4;
endRow = startRow+NRPTS_num1;
formatSpec = '
dataArray = textscan(fileID, formatSpec,NRPTS_num1+1 , 'Delimiter', delimiter,...
'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', 0,...
'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
NRPT_list = [dataArray{1:end-1}];
NRPT_list = reshape(NRPT_list',NRPTS_num1*15+15,1);
fclose(fileID);
fileID = fopen(filename,'r');
if NRPTS_num2==0
startRow = endRow;
else
startRow = endRow+1;
end
formatSpec = '
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
end
