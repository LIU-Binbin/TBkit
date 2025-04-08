function H_hr = alpharule(H_hr,level_cut,mode,options)
arguments
H_hr HR;
level_cut double{mustBeInteger} = -1;
mode double = 0;
options.Rd double = -1;
options.silence = true;
end
if options.Rd == -1
Auto_Rd = true;
else
Auto_Rd = false;
end
if mode == 0
return;
end
if level_cut == -1
level_cut = length(H_hr.Rnn);
end
Rnn = H_hr.nn_information(options.silence);
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
symvar_list_S = symvar(H_hr.ScoeL);
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
base_string_S(count_S+1:end) = [];
base_num_min_S(count_S+1:end) = [];
end
N_base_string = length(base_string);
V_n_str_L = repmat(base_string.',[level_cut-1,1]);
level_cut_str_L = string(kron((2:level_cut).',ones(N_base_string,1)));
V_n_list = sym(strcat(V_n_str_L,repmat('_',[(level_cut-1)*N_base_string,1]),level_cut_str_L),'real');
V_subs_list = sym(zeros((level_cut-1)*N_base_string,1));
if H_hr.overlap
N_base_string_S = length(base_string_S);
S_n_str_L = repmat(base_string_S.',[level_cut-1,1]);
S_n_list = sym(strcat(S_n_str_L,repmat('_',[(level_cut-1)*N_base_string_S,1]),level_cut_str_L),'real');
S_subs_list = sym(zeros((level_cut-1)*N_base_string_S,1));
end
switch mode
case 0
return;
case 1
alpha = sym('alpha');
if Auto_Rd
for n = 1:(level_cut-1)*N_base_string
[j,i] = ind2sub([level_cut-1,N_base_string],n);
j = j+1;
Coeffs_tmp =  (Rnn(j) - Rnn(base_num_min(i)))/Rnn(base_num_min(i));
V_subs_list(n) = base_symvar_min(i)*exp(-Coeffs_tmp*alpha);
end
if H_hr.overlap
alpha = sym('alpha__2');
for n = 1:(level_cut-1)*N_base_string
[j,i] = ind2sub([level_cut-1,N_base_string],n);
j = j+1;
Coeffs_tmp =  (Rnn(j) - Rnn(base_num_min_S(i)))/Rnn(base_num_min_S(i));
S_subs_list(n) = base_symvar_min_S(i)*exp(-Coeffs_tmp*alpha);
end
end
else
for n = 1:(level_cut-1)*N_base_string
[j,i] = ind2sub([level_cut-1,N_base_string],n);
j = j+1;
Coeffs_tmp =  (Rnn(j) - options.Rd)/options.Rd;
V_subs_list(n) = base_symvar_min(i)*exp(-Coeffs_tmp*alpha);
end
if H_hr.overlap
alpha = sym('alpha__2');
for n = 1:(level_cut-1)*N_base_string
[j,i] = ind2sub([level_cut-1,N_base_string],n);
j = j+1;
Coeffs_tmp =  (Rnn(j) - options.Rd)/options.Rd;
S_subs_list(n) = base_symvar_min_S(i)*exp(-Coeffs_tmp*alpha);
end
end
end
case 2
if Auto_Rd
Rd = Rnn(base_num_min);
for k = 1:length(base_symvar_min)
if isa(Rd(k),'sym')
fprintf('%s : %s \\AA\n',base_symvar_min(k),string(Rd(k)));
else
fprintf('%s : %5.3f \\AA\n',base_symvar_min(k),Rd(k));
end
end
for n = 1:(level_cut-1)*N_base_string
[j,i] = ind2sub([level_cut-1,N_base_string],n);
j = j+1;
alpha = sym(['alpha_',num2str(i),'_', char(base_symvar_min(i))]);
Coeffs_tmp =  (Rnn(j) - Rnn(base_num_min(i)))/Rnn(base_num_min(i));
V_subs_list(n) = base_symvar_min(i)*exp(-Coeffs_tmp*alpha);
end
if H_hr.overlap
Rd = Rnn(base_num_min_S);
for k = 1:length(base_symvar_min_S)
if isa(Rd(k),'sym')
fprintf('%s : %s \\AA\n',base_symvar_min_S(k),string(Rd(k)));
else
fprintf('%s : %5.3f \\AA\n',base_symvar_min_S(k),Rd(k));
end
end
for n = 1:(level_cut-1)*N_base_string
[j,i] = ind2sub([level_cut-1,N_base_string],n);
j = j+1;
alpha = sym(['alpha__2_',num2str(i),'_', char(base_symvar_min_S(i))]);
Coeffs_tmp =  (Rnn(j) - Rnn(base_num_min_S(i)))/Rnn(base_num_min_S(i));
S_subs_list(n) = base_symvar_min_S(i)*exp(-Coeffs_tmp*alpha);
end
end
else
Rd = options.Rd;
if length(Rd) == 1
Rd = repmat(Rd,[N_base_string 1]);
end
for k = 1:length(base_symvar_min)
if isa(Rd(k),'sym')
fprintf('%s : %s \\AA\n',base_symvar_min(k),string(Rd(k)));
else
fprintf('%s : %5.3f \\AA\n',base_symvar_min(k),Rd(k));
end
end
for n = 1:(level_cut-1)*N_base_string
[j,i] = ind2sub([level_cut-1,N_base_string],n);
j = j+1;
alpha = sym(['alpha_',num2str(i),'_', char(base_symvar_min(i))]);
Coeffs_tmp =  (Rnn(j) - Rd(i))/Rd(i);
V_subs_list(n) = base_symvar_min(i)*exp(-Coeffs_tmp*alpha);
end
if H_hr.overlap
for k = 1:length(base_symvar_min_S)
if isa(Rd(k),'sym')
fprintf('%s : %s \\AA\n',base_symvar_min_S(k),string(Rd(k)));
else
fprintf('%s : %5.3f \\AA\n',base_symvar_min_S(k),Rd(k));
end
end
for n = 1:(level_cut-1)*N_base_string
[j,i] = ind2sub([level_cut-1,N_base_string],n);
j = j+1;
alpha = sym(['alpha__2_',num2str(i),'_', char(base_symvar_min_S(i))]);
Coeffs_tmp =  (Rnn(j) -Rd(i))/Rd(i);
S_subs_list(n) = base_symvar_min_S(i)*exp(-Coeffs_tmp*alpha);
end
end
end
end
H_hr.HcoeL = subs(H_hr.HcoeL,V_n_list,V_subs_list);
if H_hr.overlap
H_hr.ScoeL = subs(H_hr.ScoeL,S_n_list,S_subs_list);
end
end
