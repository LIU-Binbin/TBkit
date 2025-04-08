function [simulation_result]=read_hspice_tr_sw_ac(filename)
[filepath,filetitle,ext] = fileparts(filename);
file_extension=ext(1:3);
file_handle = fopen(filename);
content_file = textscan(file_handle,'
fclose(file_handle);
content_str=string(content_file{1,1});clear content_file
data_start_ind=find(content_str==('$&
num_var=find(isnan(str2double(content_str(1:data_start_ind-2))),1)-1;
simulation_result = Hckt.save_signal_names(num_var,data_start_ind,content_str,file_extension);
data=replace(strjoin(content_str(data_start_ind:end)),' ','');
data_frmt=char(ones(1,find(char(content_str(data_start_ind))=='E',1)-1)*'.');
data_separate=double(regexp(data,[data_frmt 'E[+-]..'],'match'));clear data
for ii=1:length(simulation_result)
simulation_result(ii).val=data_separate(ii:length(simulation_result):length(data_separate)-1)';
end
end
