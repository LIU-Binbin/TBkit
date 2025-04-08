function simulation_result = save_signal_names(num_var,data_start_ind,content_str,file_extension)
expression = '[a-z]+[\d*]?+(';
TotalString = fold(@strcat,content_str);
content_char = char(TotalString);
CharIndex = regexp(content_char,expression);
InitIndex = [1,CharIndex];
EndIndex = [CharIndex-1,length(content_char)-4];
var_name_raw = repmat("",[num_var 1]);
for i = 1:num_var
var_name_raw(i) = string(content_char(InitIndex(i):EndIndex(i)));
end
var_name_raw(2:end) = var_name_raw(2:end) + ")";
if file_extension == '.tr' | file_extension == '.sw'
simulation_result=cell2struct(cellstr(var_name_raw)','var_name',1);
elseif file_extension == '.ac'
simulation_result=cell2struct(cellstr(var_name_raw)','var_name',1);
else
error(['The toolbox cannot handle' file_extension '# files']);
end
end
