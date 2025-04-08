function H_hr = full(H_hr)
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
