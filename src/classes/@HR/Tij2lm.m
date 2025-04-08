function [ml_cell,ij_list] = Tij2lm(H_hr,Rf)
NRPTS_ = H_hr.NRPTS;
orblist = H_hr.orbL;
Norb = size(orblist,1);
ij_list= unique(double(H_hr.vectorL(:,H_hr.Dim+1:H_hr.Dim+2)),'rows');
addtional_ij_list = [];
count  = 0;
ml_cell{size(ij_list,1)} = [];
for n = 1:size(ij_list,1)
t_i = orblist(ij_list(n,1),:);
t_j = orblist(ij_list(n,2),:);
t_j_m_t_i = t_j - t_i;
ml_cell{n} = [];
tmp_count = 0;
for m = 1:Norb
t_m = orblist(m,:);
for l = 1:Norb
t_l = orblist(l,:);
T_ij__ml = (t_m - t_l)*Rf - t_j_m_t_i;
if TBkit.LatticeVectorTest(T_ij__ml)
[~,j]=ismember([l,m],ij_list,'rows');
if j == 0
count = count+1;
addtional_ij_list(count,:) = [l,m];
end
tmp_count = tmp_count +1;
ml_cell{n}(tmp_count,:) = [m,l,T_ij__ml];
end
end
end
end
while ~isempty(addtional_ij_list)
ij_list = [ij_list;addtional_ij_list];
addtional_ij_list = [];
count  = 0;
for n = n+1:size(ij_list,1)
t_i = orblist(ij_list(n,1),:);
t_j = orblist(ij_list(n,2),:);
t_j_m_t_i = t_j + t_i;
ml_cell{n} = [];
tmp_count = 0;
for l = 1:Norb
t_l = orblist(l,:);
for m = 1:Norb
t_m = orblist(m,:);
T_ij__ml =  (t_m - t_l)*Rf - t_j_m_t_i;
if TBkit.LatticeVectorTest(T_ij__ml)
[~,j]=ismember([l,m],ij_list,'rows');
if j == 0
count = count+1;
addtional_ij_list(count,:) = [l,m];
end
tmp_count = tmp_count +1;
ml_cell{n}(tmp_count,:) = [l,m];
end
end
end
end
end
end
