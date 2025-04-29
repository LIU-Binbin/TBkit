function POSCAR = POSCAR_cell_read(filename,formatSpec)%POSCAR_CELL_READ Read POSCAR file into string array
%   POSCAR = POSCAR_CELL_READ(FILENAME, FORMATSPEC) reads a POSCAR file
%   and returns its contents as a string array.
%
%   Inputs:
%       filename   - Name of POSCAR file to read
%       formatSpec - Format specification for textscan
%
%   Output:
%       POSCAR     - String array containing file contents
%
%   Note:
%       Skips first row (header) and uses tab/space as delimiters
startRow = 2;
delimiter = {'\t',' '};
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
rawStringColumns = string(raw(:,:));
POSCAR=rawStringColumns;
clearvars filename delimiter startRow formatSpec fileID dataArray ans rawStringColumns raw;
end
