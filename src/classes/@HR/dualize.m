function H_hr = dualize(H_hr)
NRPTS_ = H_hr.NRPTS;
vectorList = H_hr.vectorL;
vectorList_oppo(:,1:H_hr.Dim) = -vectorList(:,1:H_hr.Dim);
if size(vectorList,2) == 5
vectorList_oppo(:,H_hr.Dim+1) = vectorList(:,H_hr.Dim+2);
vectorList_oppo(:,H_hr.Dim+2) = vectorList(:,H_hr.Dim+1);
end
for i = 1:NRPTS_
vector_tmp_oppo = vectorList_oppo(i,:);
[~,j]=ismember(vector_tmp_oppo,H_hr.vectorL,'rows');
if j == 0
H_hr = H_hr.add_empty_one(vector_tmp_oppo);
j = H_hr.NRPTS;
H_hr.Duality_vector_dist(j) = i ;
end
H_hr.Duality_vector_dist(i) = j  ;
end
end
