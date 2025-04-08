function [H_hr_R,H_hr] = applyR(H_hr,R)
arguments
H_hr HR;
R ;
end
Rf = Oper.Rc2Rf(inv(R),H_hr.Rm);
H_hr = H_hr.rewrite();
[H_hr,R_vector_dist_] = dualizeR(H_hr,Rf);
H_hr_R = H_hr;
if H_hr.num
H_hr_R.HnumL = H_hr_R.HnumL(R_vector_dist_);
end
if H_hr.coe
H_hr_R.HcoeL = H_hr_R.HcoeL(R_vector_dist_);
end
end
