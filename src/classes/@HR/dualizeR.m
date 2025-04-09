function [H_hr,R_vector_dist_] = dualizeR(H_hr,Rf)
% DUALIZER Apply dual R operation to Hamiltonian in HR format
%
%   [H_hr,R_vector_dist_] = DUALIZER(H_hr,Rf) transforms the Hamiltonian
%   using the given R matrix and returns the modified Hamiltonian along
%   with distance vectors.
%
%   INPUT ARGUMENTS:
%       H_hr - Hamiltonian in HR format (HR object)
%       Rf - Transformation matrix (3x3 double)
%
%   OUTPUT ARGUMENTS:
%       H_hr - Modified Hamiltonian in HR format
%       R_vector_dist_ - Distance vectors between transformed points
%
%   NOTES:
%       - Modifies the input HR object by adding empty points when needed
%       - Uses floor operation for vector calculations
%
%   SEE ALSO:
%       HR, timtj_gen
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]
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
