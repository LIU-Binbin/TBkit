function [H_hr,R_vector_dist_] = dualizeR(H_hr,Rf)
NRPTS_ = H_hr.NRPTS;
ONESLIST = ones(NRPTS_,1);
H_hr = H_hr.timtj_gen();
timtj_mat = H_hr.timtj{2};
Size_timtj_mat = size(timtj_mat);
vectorList = double(H_hr.vectorL(:,1:H_hr.Dim));
IndList1 = sub2ind(Size_timtj_mat,H_hr.vectorL(:,H_hr.Dim+1),H_hr.vectorL(:,H_hr.Dim+2),ONESLIST);
IndList2 = sub2ind(Size_timtj_mat,H_hr.vectorL(:,H_hr.Dim+1),H_hr.vectorL(:,H_hr.Dim+2),ONESLIST*2);
IndList3 = sub2ind(Size_timtj_mat,H_hr.vectorL(:,H_hr.Dim+1),H_hr.vectorL(:,H_hr.Dim+2),ONESLIST*3);
vectorL_addtional = [timtj_mat(IndList1),timtj_mat(IndList2),timtj_mat(IndList3)];
vectorL_  = floor(H_hr.orbL(H_hr.vectorL(:,H_hr.Dim+1),:)+(vectorList-vectorL_addtional)*double(Rf));
for i = 1:NRPTS_
vector_tmp_oppo = vectorL_(i,:);
[~,j]=ismember(vector_tmp_oppo,vectorList,'rows');
if j == 0
H_hr = H_hr.add_empty_one([vector_tmp_oppo,H_hr.vectorL(i,14:5)]);
j = H_hr.NRPTS;
H_hr.R_vector_dist(j) = i ;
end
H_hr.R_vector_dist(i) = j  ;
R_vector_dist_ = H_hr.R_vector_dist;
end
end
