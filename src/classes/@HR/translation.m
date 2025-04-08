function H_hr = translation(H_hr,translation_vector,options)
arguments
H_hr HR;
translation_vector double = [0,0,0];
options.Accuracy double = 1e-6;
options.force_list = false;
end
OUT_H_hr = H_hr;
pc_orb = H_hr.orbL;
OUT_H_hr.orbL = mod(H_hr.orbL + translation_vector,1);
sc_orb = OUT_H_hr.orbL;
if strcmp(H_hr.Type,'sqarse')
H_hr = H_hr.full();
end
if strcmp(H_hr.Type,'mat')
fprintf("Translation must use list mode!");
H_hr = H_hr.rewrite();
OUT_H_hr = OUT_H_hr.rewrite();
GiveBackMat =true;
else
GiveBackMat =false;
end
if H_hr.num
HnumList = H_hr.HnumL;
end
if H_hr.coe
HcoeList = H_hr.HcoeL;
end
VectorList = double(H_hr.vectorL);
hiL = VectorList(:,H_hr.Dim+1);
hjL = VectorList(:,H_hr.Dim+2);
ind_RL = double(VectorList(:,1:H_hr.Dim));
sc_hjL = hjL;
sc_hiL = hiL;
indRtiL = pc_orb(hiL,:);
indRti_in_supercellL = sc_orb(sc_hiL,:);
real_sc_vecL = indRti_in_supercellL - indRtiL;
indRtjL = real_sc_vecL+ind_RL+pc_orb(hjL,:);
indRtj_in_supercellL = (double(indRtjL));
indR_in_supercellL  = floor(indRtj_in_supercellL);
OutVectorList= ...
[indR_in_supercellL,sc_hiL,sc_hjL];
if H_hr.num
OUT_H_hr.HnumL = HnumList;
end
if H_hr.coe
OUT_H_hr.HcoeL = HcoeList;
end
OUT_H_hr.vectorL = OutVectorList;
H_hr = OUT_H_hr;
H_hr.Basis_num = OUT_H_hr.WAN_NUM;
if options.force_list
if strcmp(H_hr.Type,'mat')
H_hr = H_hr.rewrite();
end
elseif GiveBackMat
if ~strcmp(H_hr.Type,'mat')
H_hr = H_hr.rewind();
end
else
end
end
