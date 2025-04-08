function [H_htrig,H_htrig2] = applyR(H_htrig,R)
arguments
H_htrig Htrig;
R ;
end
[num_label,coe_label] = H_htrig.NumOrCoe();
[H_htrig,Smat] = H_htrig.Smatgen(R);
H_htrig2 = H_htrig;
if num_label
H_htrig.HnumL = Htrig.matrixtimespage(Smat,H_htrig.HnumL);
end
if coe_label
H_htrig.HcoeL = Htrig.matrixtimespage(Smat,H_htrig.HcoeL);
end
end
