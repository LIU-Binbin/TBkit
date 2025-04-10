function H_hr = nn_sk_smart(H_hr,search_range,Accuracy,Rlength_cut)
% NN_SK_SMART Generate smart Slater-Koster neighbor information
%
%   H_HR = NN_SK_SMART(H_HR,SEARCH_RANGE,ACCURACY,RLENGTH_CUT) generates
%   smart Slater-Koster neighbor tables for HR object
%
%   Inputs:
%       H_hr - HR object to populate with NN information
%       search_range - Search range vector [default: [0 0 0]]
%       Accuracy - Rounding accuracy [default: 4]
%       Rlength_cut - Cutoff distance [default: 15]
%   Output:
%       H_hr - Updated HR object with neighbor information
%
%   Notes:
%       - Builds both Atom_store_smart and nn_store_smart structures
%       - Creates Rnn distance list and mapping container
%       - Includes progress display during generation
if nargin <4
Rlength_cut = 15;
end
if nargin <3
Accuracy = 4;
end
if nargin <2
search_range = [0 0 0];
end
sites_num = size(H_hr.sites,2);
search_rangex=search_range(1);
search_rangey=search_range(2);
search_rangez=search_range(3);
Atom_smart_t = struct('R_fractional_from',[],...
'R_fractional_to',[],'R_fractional_diff',[],...
'seq_from',[],'seq_to',[],...
'orb_sym_from',[],'orb_sym_to',[],...
'l_name_from',[],'l_name_to',[],...
'handyname',[]);
nn_smart_t = struct('seq_from',[],'seq_to',[],'nn',[]);
H_hr.Atom_store_smart = repmat(Atom_smart_t,[sites_num sites_num]);
H_hr.nn_store_smart = repmat(nn_smart_t,[sites_num sites_num]);
Rnn_list = [];
pb = TBkit_tool_outer.CmdLineProgressBar('Seaching ');
for j  = 1:sites_num
site2 = H_hr.sites(j);
for i = 1:sites_num
site1 = H_hr.sites(i);
H_hr.Atom_store_smart(i,j) = HR.Atom_smart_t_gen(site1,site2);
[Rnn_list_temp,H_hr.nn_store_smart(i,j)] = ...
HR.nn_smart_t_gen(H_hr.Atom_store_smart(i,j),H_hr.Rm,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut);
Rnn_list = [Rnn_list;Rnn_list_temp];
end
pb.print(j,sites_num,' th orb nn information ...');
end
pb.delete();
H_hr.Rnn = sort(unique(Rnn_list,'row'));
H_hr.Rnn_map = containers.Map('KeyType','double','ValueType','double');
for i = 1:length(H_hr.Rnn)
H_hr.Rnn_map(H_hr.Rnn(i)) = i;
end
pb = TBkit_tool_outer.CmdLineProgressBar('Giving ');
for j  = 1:sites_num
pb.print(j,sites_num,' th orb nn_level ...');
for i = 1:sites_num
for k = 1:length(H_hr.nn_store_smart(i,j).nn)
nn_level = H_hr.Rnn_map(H_hr.nn_store_smart(i,j).nn(k).Rlength);
H_hr.nn_store_smart(i,j).nn(k).nn_level = nn_level ;
H_hr.nn_store_smart(i,j).nn(k).hop = HR.hop_gen(H_hr.nn_store_smart(i,j).nn(k).hop_pre,nn_level);
if H_hr.overlap
H_hr.nn_store_smart(i,j).nn(k).overlap = sym(strrep(string(H_hr.nn_store_smart(i,j).nn(k).hop),'V','S'));
end
end
end
end
pb.delete();
end
