function H_hr = GenfromOrth(H_hr,seed_r,seed_i,Accuracy,options)
arguments
H_hr HR;
seed_r = 'gamma__r_';
seed_i = 'gamma__i_';
Accuracy = 1e-6;
options.fromCvectorL = false;
end
nAccuracy = floor(-log(Accuracy)/log(10));
switch H_hr.Type
case 'list'
if H_hr.vectorhopping
if  options.fromCvectorL
CL= vpa(H_hr.CvectorL(:,1:rank(H_hr.CvectorL)),nAccuracy);
SymVar_r = sym(seed_r,[size(CL,2),1],'real');
H_hr.HcoeL = CL(1:end/2,:)*SymVar_r + 1i*CL(end/2+1:end,:)*SymVar_r;
H_hr.vectorhopping = false;
return;
end
AL = vpa(H_hr.AvectorL,nAccuracy);
BL = vpa(H_hr.BvectorL,nAccuracy);
SymVar_r = sym(seed_r,[size(AL,2),1],'real');
SymVar_i = sym(seed_i,[size(BL,2),1],'real');
H_hr.HcoeL = AL*SymVar_r + 1i*BL*SymVar_i;
H_hr.vectorhopping = false;
end
end
end
