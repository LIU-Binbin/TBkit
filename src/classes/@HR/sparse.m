function H_hr = sparse(H_hr)
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
