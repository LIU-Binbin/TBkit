function [H_hr_R,H_hr] = applyRU(H_hr,SymOper )
arguments
H_hr HR;
SymOper     ;
end
[H_hr,VectorDistMat] = dualizeOper(H_hr,SymOper);
if size(VectorDistMat,1)~=size(VectorDistMat,2)
error('check why?');
end
H_hr_R = H_hr;
if H_hr.vectorhopping
AvectorLtmp = H_hr_R.AvectorL;
BvectorLtmp = H_hr_R.BvectorL;
CvectorLtmp = H_hr_R.CvectorL;
if SymOper.conjugate
BvectorLtmp = -(BvectorLtmp);
CvectorLtmp(end/2+1:end,:) = -CvectorLtmp(end/2+1:end,:);
end
if SymOper.antisymmetry
AvectorLtmp = -AvectorLtmp;
BvectorLtmp = -BvectorLtmp;
CvectorLtmp = -CvectorLtmp;
end
H_hr_R.AvectorL = VectorDistMat*AvectorLtmp;
H_hr_R.BvectorL = VectorDistMat*BvectorLtmp;
CL1 = VectorDistMat*CvectorLtmp(1:end/2,:);
CL2 = VectorDistMat*CvectorLtmp(end/2+1:end,:);
H_hr_R.CvectorL = [real(CL1)-imag(CL2);imag(CL1)+real(CL2)];
return;
end
if H_hr.coe
HcoeLtmp = H_hr_R.HcoeL;
if SymOper.conjugate
HcoeLtmp = conj(HcoeLtmp);
end
if SymOper.antisymmetry
HcoeLtmp = -HcoeLtmp;
end
H_hr_R.HcoeL = VectorDistMat*HcoeLtmp;
end
if H_hr.num == true
HnumLtmp = H_hr_R.HnumL;
if SymOper.conjugate
HnumLtmp = conj(HnumLtmp);
end
if SymOper.antisymmetry
HnumLtmp = -HnumLtmp;
end
H_hr_R.HnumL = VectorDistMat*HnumLtmp;
end
end
