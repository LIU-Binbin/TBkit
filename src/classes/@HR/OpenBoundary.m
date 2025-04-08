function H_hr = OpenBoundary(H_hr,OBC_list)
arguments
H_hr ;
OBC_list = zeros(1,H_hr.Dim);
end
if isequal(OBC_list,zeros(1,H_hr.Dim))
return;
end
Type = H_hr.Type;
H_hr =H_hr.ForceToMat();
vectorList = H_hr.vectorL;
ABSvectorList = abs(vectorList);
KeepVectorL = find(~logical(sum(ABSvectorList(:,logical(OBC_list))>0,2)));
H_hr = H_hr.reseq(':',KeepVectorL);
H_hr = H_hr.ForceToType(Type);
end
