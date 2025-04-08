function H_hr = project(H_hr,BASIS_MAT)
new_WAN_NUM = size(BASIS_MAT,1);
BASIS_MAT_prime = BASIS_MAT';
switch H_hr.Type
case 'sparse'
disp('waiting');
case 'list'
disp('waiting');
case 'mat'
if H_hr.num
H_hr.HcoeL = zeros(new_WAN_NUM,new_WAN_NUM,H_hr.NRPTS);
tempHnumL = zeros(new_WAN_NUM,new_WAN_NUM,H_hr.NRPTS);
for i = 1:H_hr.NRPTS
tempHnumL(:,:,i) = BASIS_MAT*H_hr.HnumL(:,:,i)*BASIS_MAT_prime ;
end
H_hr.HnumL = tempHnumL;
end
if H_hr.coe
H_hr.HnumL = zeros(new_WAN_NUM,new_WAN_NUM,H_hr.NRPTS);
tempHcoeL = sym(zeros(new_WAN_NUM,new_WAN_NUM,H_hr.NRPTS));
for i = 1:H_hr.NRPTS
tempHcoeL(:,:,i) = BASIS_MAT*H_hr.HcoeL(:,:,i)*BASIS_MAT_prime ;
end
H_hr.HcoeL = tempHcoeL;
end
end
end
