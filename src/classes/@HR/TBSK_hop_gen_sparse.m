function [TBSK_hop,Coff] =  TBSK_hop_gen_sparse(site1,site2,Rij_cart,Rlength,nn_level)
Rlmn = Rij_cart/Rlength;
orb1 = site1.orb   ;
orb2 = site2.orb   ;
orb_sym1 =  site1.orb_sym;
orb_sym2 =  site2.orb_sym;
orbsym1_n =  HR.subs_xyz(orb_sym1 ,Rlmn);
orbsym2_n =  HR.subs_xyz(orb_sym2 ,Rlmn);
if strcmp(orb1,'p')
if strcmp(orb2,'s')
orb1 = 's';
orb2 = 'p';
orbsym1_n = -orbsym1_n;
end
end
Coff(1) = orbsym1_n*orbsym2_n;
Coff(2) = HR.delta_orb(orb1,orb2)*(HR.delta_orb_sym(orb_sym1,orb_sym2)-orbsym1_n*orbsym2_n);
if abs(Coff(1)) < 1e-10
Coff(1) = 0;
end
if abs(Coff(2)) < 1e-10
Coff(2) = 0;
end
Coff(3) = 0;
VP = "V"+orb1+orb2+'S_'+num2str(nn_level);
VS = "V"+orb1+orb2+'P_'+num2str(nn_level);
VD = "";
TBSK_hop =  Coff(1) * sym(VP,'real')+...
Coff(2) * sym(VS,'real');
end
