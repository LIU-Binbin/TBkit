function H_hr = simplify(H_hr,Accuracy)
% SIMPLIFY Simplify HR object by removing negligible terms
%
%   H_HR = SIMPLIFY(H_HR,ACCURACY) simplifies HR object contents
%
%   Inputs:
%       H_hr - HR object to simplify
%       Accuracy - Threshold for term removal [default: 1e-6]
%   Output:
%       H_hr - Simplified HR object
%
%   Notes:
%       - Removes terms below accuracy threshold
%       - Handles vector hopping terms specially
%       - Works on both numeric and symbolic coefficients
%       - Maintains storage format consistency
if nargin < 2
Accuracy = 1e-6;
end
if H_hr.vectorhopping
AL = rref(H_hr.AvectorL.').';
BL = rref(H_hr.BvectorL.').';
CL = rref(H_hr.CvectorL.').';
AL = AL(:,1:rank(AL));
BL = BL(:,1:rank(BL));
CL = CL(:,1:rank(CL));
rAL = real(AL);rAL(abs(rAL)<Accuracy) = 0;
rBL = real(BL);rBL(abs(rBL)<Accuracy) = 0;
rCL = real(CL);rCL(abs(rCL)<Accuracy) = 0;
H_hr.AvectorL = rAL;
H_hr.BvectorL = rBL;
H_hr.CvectorL = rCL;
return;
end
if H_hr.coe
H_coeL_tmp = simplify(H_hr.HcoeL);
H_hr.HcoeL = H_coeL_tmp;
if strcmp(H_hr.Type,'list')
NRPTS_list = find(H_coeL_tmp ~=sym(0));
H_hr = H_hr.reseq(':',NRPTS_list);
elseif strcmp(H_hr.Type,'mat')
end
end
if H_hr.num
H_numL_tmp = H_hr.HnumL;
if strcmp(H_hr.Type,'list')
NRPTS_list = find(abs(H_numL_tmp) > Accuracy);
H_hr = H_hr.reseq(':',NRPTS_list);
elseif strcmp(H_hr.Type,'mat')
zerosMat = ones(size(H_numL_tmp(:,:,1)))*Accuracy;
NRPTS_list = true(H_hr.NRPTS,1);
for i = 1:H_hr.NRPTS
if sum(abs(H_numL_tmp(:,:,i)) > zerosMat,'all')
else
NRPTS_list(i) = false;
end
end
H_hr = H_hr.reseq(':',NRPTS_list);
end
end
end
