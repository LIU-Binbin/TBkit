function [H_hr_forHckt,maxOnsite] = HRforHckt(H_hr,options)
arguments
H_hr HR;
options.fast = false
options.direct = false;
options.Accuracy = 1e-6;
options.C_0 = 1;
options.coefficient = -1;
end
if strcmp(H_hr.Type,'mat')
H_hr = H_hr.rewrite();
end
H_hr_forHckt = H_hr;
if H_hr.coe
H_hr_forHckt.HcoeL = options.coefficient * H_hr_forHckt.HcoeL;
C_0 = sym('C_0','real');
assume(C_0>0);
BaseOnsiteL = repmat(C_0,[H_hr_forHckt.WAN_NUM,1]);
maxOnsite = C_0;
for i = 1:H_hr_forHckt.WAN_NUM
maxOnsite = max(BaseOnsiteL(i)-sum(H_hr_forHckt.HcoeL(H_hr_forHckt.vectorL(:,H_hr.Dim+1)==i)),maxOnsite);
end
%fprintf('The universal shift is force set to be: %c',maxOnsite);
for i = 1:H_hr_forHckt.WAN_NUM
BaseOnsiteL(i) = maxOnsite + sum(H_hr_forHckt.HcoeL(H_hr_forHckt.vectorL(:,H_hr.Dim+1)==i));
H_hr_forHckt = H_hr_forHckt.set_hop(maxOnsite,...
i,i,zeros(1,H_hr.Dim),'symadd');
end
else
H_hr_forHckt.HnumL = options.coefficient * H_hr_forHckt.HnumL;
C_0 = options.C_0 ;
BaseOnsiteL = repmat(C_0,[H_hr_forHckt.WAN_NUM,1]);
maxOnsite = C_0;
for i = 1:H_hr_forHckt.WAN_NUM
maxOnsite = max(BaseOnsiteL(i)-sum(H_hr_forHckt.HnumL(H_hr_forHckt.vectorL(:,H_hr.Dim+1)==i)),maxOnsite);
end
% fprintf('The universal shift is forcely set to be: %f',maxOnsite);
for i = 1:H_hr_forHckt.WAN_NUM
BaseOnsiteL(i) = maxOnsite + sum(H_hr_forHckt.HnumL(H_hr_forHckt.vectorL(:,H_hr.Dim+1)==i));
H_hr_forHckt = H_hr_forHckt.set_hop(maxOnsite,...
i,i,zeros(1,H_hr.Dim),'add');
end
end
end
