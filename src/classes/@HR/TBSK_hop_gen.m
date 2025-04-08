function [TBSK_hop,Coff] = TBSK_hop_gen(orb1,orb2,orbsym1_n,orbsym2_n,orbsym1,orbsym2)
if strcmp(orb1,'p')
if strcmp(orb2,'s')
orb1 = 's';
orb2 = 'p';
orbsym1_n = -orbsym1_n;
end
end
Coff(1) = orbsym1_n*orbsym2_n;
Coff(2) = HR.delta_orb(orb1,orb2)*(HR.delta_orb_sym(orbsym1,orbsym2)-orbsym1_n*orbsym2_n);
Coff(3) = 0;
TBSK_hop = Coff(1) *str2sym("V"+orb1+orb2+'S')+...
Coff(2)*str2sym("V"+orb1+orb2+'P');
end
