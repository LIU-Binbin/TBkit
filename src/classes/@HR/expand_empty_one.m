function H_hr = expand_empty_one(H_hr,orbOne,QuantumOne,elementOne)
arguments
H_hr;
orbOne = zeros(1,H_hr.Dim);
QuantumOne = [1,0,0,1];
elementOne = 1;
end
H_hr.orbL = [H_hr.orbL;orbOne];
H_hr.quantumL = [H_hr.quantumL;QuantumOne];
H_hr.elementL = [H_hr.elementL;elementOne];
if strcmp(H_hr.Type ,'list')
NRPTS_new = H_hr.NRPTS +size(orbOne,1);
H_hr.HcoeL(NRPTS_new,1) = sym(0);
H_hr.HnumL(NRPTS_new,1) = 0;
else
WANNUM = H_hr.WAN_NUM+size(orbOne,1);
H_hr.HcoeL(WANNUM,WANNUM,:) = sym(0);
H_hr.HnumL(WANNUM,WANNUM,:) = 0;
end
end
