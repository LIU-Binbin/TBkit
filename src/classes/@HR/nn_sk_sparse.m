function H_hr = nn_sk_sparse(H_hr,search_range,Accuracy,Rlength_cut)
% NN_SK_SPARSE Generate sparse Slater-Koster neighbor information
%
%   H_HR = NN_SK_SPARSE(H_HR,SEARCH_RANGE,ACCURACY,RLENGTH_CUT) generates
%   sparse Slater-Koster neighbor tables for HR object
%
%   Inputs:
%       H_hr - HR object to populate with NN information
%       search_range - Search range vector [default: [0 0 0]]
%       Accuracy - Rounding accuracy [default: 4]
%       Rlength_cut - Cutoff distance [default: 5]
%   Output:
%       H_hr - Updated HR object with neighbor information
%
%   Notes:
%       - Uses sparse matrix storage format
%       - More memory-efficient for large systems
%       - Includes progress display during generation
MAX_NN_LENGTH = 5000*5000*9;
if nargin <4
Rlength_cut = 5;
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
max_nn_length = min(H_hr.WAN_NUM*H_hr.WAN_NUM*H_hr.NRPTS,MAX_NN_LENGTH);
nn_sparse = zeros(max_nn_length,10);
Rnn_list = [];
count  = 1;
for j  = 1:sites_num
site2 = H_hr.sites(j);
for i = 1:sites_num
site1 = H_hr.sites(i);
[nn_sparse_temp,Rnn_list_temp] = ...
HR.nn_sparse_t_gen(site1,site2,H_hr.Rm,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut);
countadd = size(nn_sparse_temp,1);
nn_sparse_temp(:,1) = i;
nn_sparse_temp(:,2) = j;
nn_sparse(count:count+countadd-1,:) = nn_sparse_temp;
count = count+countadd;
Rnn_list = [Rnn_list;Rnn_list_temp];
end
fprintf('The (%4d/%4d) th orb nn information has been generated.\n',j,sites_num);
end
nn_sparse(count:max_nn_length,:) = [];
H_hr.Rnn = sort(unique(Rnn_list,'row'));
H_hr.Rnn_map = containers.Map('KeyType','double','ValueType','double');
for i = 1:length(H_hr.Rnn)
H_hr.Rnn_map(H_hr.Rnn(i)) = i;
end
for i  = 1:count-1
fprintf('Giving (%4d/%4d) th hopping  ... \n',i,count);
nn_sparse(i,7) = H_hr.Rnn_map(nn_sparse(i,6));
end
H_hr.nn_store_smart = nn_sparse;
end
