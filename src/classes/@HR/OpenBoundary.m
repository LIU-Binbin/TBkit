function H_hr = OpenBoundary(H_hr,OBC_list)
% OPENBOUNDARY Apply open boundary conditions to HR object
%
%   H_HR = OPENBOUNDARY(H_HR,OBC_LIST) applies open boundary conditions
%   along specified directions
%
%   Inputs:
%       H_hr - HR object to modify
%       OBC_list - List of directions to open [default: zeros]
%   Output:
%       H_hr - Modified HR object with open boundaries
%
%   Notes:
%       - Preserves original matrix type
%       - No effect if OBC_list is all zeros
%       - Forces to matrix form temporarily during processing
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
