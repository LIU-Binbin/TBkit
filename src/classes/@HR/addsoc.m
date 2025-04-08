function H_hr = addsoc(H_hr,quantumL)
if nargin >1
H_hr.quantumL = quantumL;
end
H_hr = H_hr + H_hr.H_atom_soc;
end
