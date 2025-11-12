function DENSITY_SET = CHGCAR_raw_read(filename)
arguments
    filename
end
%% find the FFT grid = NGXF*NGYF*NGZF 
fid = fopen(filename, 'r');
while true
    fline = fgetl(fid);
    if isempty(strtrim(fline))
        break;
    end
end
FFT_grid_line = fgetl(fid);
tmp_list = split(string(FFT_grid_line));
NGXF = str2double(tmp_list(2));
NGYF = str2double(tmp_list(3));
NGZF = str2double(tmp_list(4));
%% find the start of each data block
frewind(fid);
line_number = 0;
matches = [];
while ~feof(fid)
    line_number = line_number +1;
    fline = fgetl(fid);
    if strcmp(fline, FFT_grid_line)
        matches = [matches, line_number];
    end
end
starts = matches +1;
fclose(fid);
%% read each data block
NGFtot = NGXF * NGYF * NGZF;
block_length = ceil(NGFtot/5);

DENSITY_SET = cell(length(starts),1);
for ib = 1:length(starts)
    fid = fopen(filename,'r');
    startline = starts(ib);
    formatSpec = '%f%f%f%f%f';
    datablock = textscan(fid, formatSpec, block_length, 'Delimiter', ' ',...
        'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', 0,...
        'HeaderLines', startline-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    density = [datablock{1:end}];
    density = reshape(density', [], 1);
    density = reshape(density(1:NGFtot), NGXF, NGYF, NGZF);
    fclose(fid);

    DENSITY_SET{ib} = density;
end
end
