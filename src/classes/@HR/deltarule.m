function H_hr = deltarule(H_hr,level_cut,mode,options)
arguments
H_hr HR;
level_cut double{mustBeInteger} = 1;
mode double = 0;
options.Rd = -1;
end
import park.*;
if mode == 0
return;
end
if level_cut == -1
level_cut = length(H_hr.Rnn);
end
Rnn = H_hr.nn_information();
base_string = ["VssS","VspS","VsdS","VppS","VpdS","VppP","VpdP","VddS","VddP","VddD"];
strvar_list = string(H_hr.symvar_list);
base_symvar_min = sym([]);
base_num_min = ones(1,length(base_string));
count = 0;
for ibase_string = base_string
for i = 1:level_cut
if contains(ibase_string+"_"+string(i),strvar_list)
count = count +1;
base_symvar_min(count) = sym(ibase_string+"_"+string(i),'real');
base_num_min(count) = i;
base_string(count) = ibase_string;
break;
end
end
end
base_string(count+1:end) = [];
base_num_min(count+1:end) = [];
if H_hr.overlap
base_string_S = ["SssS","SspS","SsdS","SppS","SpdS","SppP","SpdP","SddS","SddP","SddD"];
symvar_list_S = symvar(TBkitobj.ScoeL);
strvar_list_S  = [string(symvar_list_S)];
base_symvar_min_S = sym([]);
base_num_min_S = ones(1,length(base_string));
count_S = 0;
for ibase_string_S = base_string_S
for i = 1:level_cut
if contains(ibase_string_S+"_"+string(i),strvar_list_S)
count_S = count_S +1;
base_symvar_min_S(count_S) = sym(ibase_string_S+"_"+string(i),'real');
base_num_min_S(count_S) = i;
base_string_S(count_S) = ibase_string_S;
break;
end
end
end
base_string_S(count+1:end) = [];
base_num_min_S(count+1:end) = [];
end
if options.Rd == -1
RdL = Rnn;
else
disp(base_string);
if length(options.Rd) > 1
RdL = repmat(options.Rd,[1 count]);
else
RdL = options.Rd;
end
base_num_min = 1:count;
if H_hr.overlap
base_num_min_S = 1:count_S;
end
end
switch mode
case 0
return;
case 1
for j = 2:level_cut
for i= 1:length(base_string)
ibase_string = base_string(i);
V_n = ibase_string+"_"+string(j);
if strcontain(V_n ,strvar_list)
delta = sym('delta','real');
Coeffs_tmp =  (Rnn(j) - RdL(base_num_min(i)));
V_subs = base_symvar_min(i)*exp(-Coeffs_tmp/delta);
H_hr.HcoeL = subs(H_hr.HcoeL,sym(V_n),V_subs);
if H_hr.overlap
H_hr.ScoeL = subs(H_hr.ScoeL,sym(V_n),V_subs);
end
end
end
end
if H_hr.overlap
for j = 2:level_cut
for i= 1:length(base_string_S)
ibase_string_S = base_string_S(i);
S_n = ibase_string_S+"_"+string(j);
if strcontain(S_n ,strvar_list_S)
delta = sym('delta__2','real');
Coeffs_tmp =  (Rnn(j) - RdL(base_num_min_S(i)));
S_subs = base_num_min_S(i)*exp(-Coeffs_tmp/delta);
H_hr.ScoeL = subs(H_hr.ScoeL,sym(S_n),S_subs);
end
end
end
end
case 2
for j = 2:level_cut
for i= 1:length(base_string)
ibase_string = base_string(i);
V_n = ibase_string+"_"+string(j);
if strcontain(V_n ,strvar_list)
delta = sym(['delta_',num2str(i)],'real');
Coeffs_tmp =  (Rnn(j) - RdL(base_num_min(i)));
V_subs = base_symvar_min(i)*exp(-Coeffs_tmp/delta);
H_hr.HcoeL = subs(H_hr.HcoeL,sym(V_n),V_subs);
end
end
end
if H_hr.overlap
for j = 2:level_cut
for i= 1:length(base_string_S)
ibase_string_S = base_string_S(i);
S_n = ibase_string_S+"_"+string(j);
if strcontain(S_n ,strvar_list_S)
delta = sym(['delta__2_',num2str(i)],'real');
Coeffs_tmp =  (Rnn(j) - RdL(base_num_min_S(i)));
S_subs = base_num_min_S(i)*exp(-Coeffs_tmp/delta);
H_hr.ScoeL = subs(H_hr.ScoeL,sym(S_n),S_subs);
end
end
end
end
end
end
