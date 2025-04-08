function H_hr = add_soc(H_hr)
if ~H_hr.soc
quantumList_up = H_hr.quantumL;
quantumList_up(:,4) = 1;
quantumList_dn = quantumList_up;
quantumList_dn(:,4) = -1;
H_hr = [H_hr,(H_hr)];
end
H_hr.quantumL = [quantumList_up;quantumList_dn];
H_hr.HcoeL = sym(H_hr.HnumL);
H_hr.coe = true;
H_hr.num = false;
H_hr = H_hr + H_hr.H_atom_soc();
end
