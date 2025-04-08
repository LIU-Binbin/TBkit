function [H_hr_out,H_hr_pi_plus,H_hr_pi_minus] = realmap(H_hr)
switch H_hr.Type
case 'sparse'
H_hr = H_hr.full;
case 'list'
H_hr = H_hr.rewind;
case 'mat'
end
H_hr_out = H_hr;
H_hr_out.quantumL = [H_hr_out.quantumL;H_hr_out.quantumL];
H_hr_out.elementL = [H_hr_out.elementL;H_hr_out.elementL];
H_hr_out.orbL = [H_hr_out.orbL;H_hr_out.orbL];
new_WAN_NUM = 2*H_hr_out.WAN_NUM;
WANNUM = H_hr.WAN_NUM;
NRPTS_ = H_hr.NRPTS;
PI_plus = 1/2 * [eye(WANNUM), -1i*eye(WANNUM);1i*eye(WANNUM),eye(WANNUM)];
PI_minus = 1/2 * [eye(WANNUM), 1i*eye(WANNUM);-1i*eye(WANNUM),eye(WANNUM)];
if H_hr.coe
HcoeList = H_hr.HcoeL;
realHcoeList = real(HcoeList);
imagHcoeList = imag(HcoeList);
tempHcoeL = sym(zeros(new_WAN_NUM,new_WAN_NUM,H_hr.NRPTS));
tempHcoeL_pi_plus = tempHcoeL;
tempHcoeL_pi_minus = tempHcoeL;
for i = 1:NRPTS_
tempHcoeL(:,:,i) = [realHcoeList(:,:,i),imagHcoeList(:,:,i);-imagHcoeList(:,:,i),realHcoeList(:,:,i)] ;
tempHcoeL_pi_plus(:,:,i)  = PI_plus*tempHcoeL(:,:,i)*PI_plus;
tempHcoeL_pi_minus(:,:,i)  = PI_minus*tempHcoeL(:,:,i)*PI_minus;
end
H_hr_out.HcoeL = tempHcoeL;
H_hr_pi_plus = H_hr_out;
H_hr_pi_plus.HcoeL = tempHcoeL_pi_plus;
H_hr_pi_minus = H_hr_out;
H_hr_pi_minus.HcoeL = tempHcoeL_pi_minus;
elseif H_hr.num
HnumList = H_hr.HnumL;
realHnumList = real(HnumList);
imagHnumList = imag(HnumList);
tempHnumL = sym(zeros(new_WAN_NUM,new_WAN_NUM,H_hr.NRPTS));
PI_plus_page = repmat(PI_plus,[1 1 NRPTS_]);
PI_minus_page = repmat(PI_minus,[1 1 NRPTS_]);
for i = 1:HNRPTS_
tempHnumL(:,:,i) = [realHnumList(:,:,i),imagHnumList(:,:,i);-imagHnumList(:,:,i),realHnumList(:,:,i)] ;
end
H_hr_out.HnumL = tempHnumL;
H_hr_pi_plus = H_hr_out;
H_hr_pi_plus.HnumL = pagemtimes(pagemtimes(PI_plus_page,tempHnumL),PI_plus_page);
H_hr_pi_minus = H_hr_out;
H_hr_pi_minus.HnumL = pagemtimes(pagemtimes(PI_minus_page,tempHnumL),PI_minus_page);
end
end
