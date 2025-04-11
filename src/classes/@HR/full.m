function H_hr = full(H_hr)
% FULL Convert sparse HR object to dense matrix format
%
%   H_hr = FULL(H_hr) converts a sparse HR object to dense matrix format.
%
%   INPUT ARGUMENTS:
%       H_hr - Sparse HR object to convert
%
%   OUTPUT ARGUMENTS:
%       H_hr - Dense matrix format HR object
%
%   NOTES:
%       - Only converts if Type is 'sparse'
%       - Returns input unchanged for other types
%
%   SEE ALSO:
%       HR, sparse
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

if ~strcmp(H_hr.Type,'sparse')
return;
end
HnumL_temp = zeros(H_hr.WAN_NUM,H_hr.WAN_NUM,H_hr.NRPTS);
for i = 1:H_hr.NRPTS
HnumL_temp(:,:,i) = full(H_hr.HnumL{i});
end
H_hr.HnumL = HnumL_temp;
H_hr.Type = 'mat';
end
