function H_hr = sparse(H_hr)
% SPARSE Convert HR object to sparse storage format
%
%   H_HR = SPARSE(H_HR) converts HR object to sparse format
%
%   Input:
%       H_hr - HR object to convert
%   Output:
%       H_hr - HR object in sparse format
%
%   Notes:
%       - Only converts numeric terms
%       - Clears symbolic coefficients
%       - Maintains all hopping information
%       - No effect if already in sparse format
if strcmp(H_hr.Type,'sparse')
return;
end
HnumL_temp{H_hr.NRPTS} = zeros(H_hr.WAN_NUM,H_hr.WAN_NUM);
for i = 1:H_hr.NRPTS
HnumL_temp{i}= sparse(H_hr.HnumL(:,:,i));
end
H_hr.HnumL = HnumL_temp;
H_hr.HcoeL = [];
H_hr.Type = 'sparse';
end
