function [TBSK_hop,Coff] = TBSK_hop_gen(orb1,orb2,orbsym1_n,orbsym2_n,orbsym1,orbsym2)
%TBSK_HOP_GEN Generate Slater-Koster hopping terms
%
%   [TBSK_HOP, COFF] = TBSK_HOP_GEN(ORB1, ORB2, ORBSYM1_N, ORBSYM2_N, ORBSYM1, ORBSYM2)
%   Generates Slater-Koster hopping terms and coefficients for tight-binding calculations
%
%   Inputs:
%       orb1      - Orbital type of first atom ('s', 'p', 'd', 'f')
%       orb2      - Orbital type of second atom ('s', 'p', 'd', 'f')
%       orbsym1_n - Direction cosine component for first orbital
%       orbsym2_n - Direction cosine component for second orbital
%       orbsym1   - Symbolic representation of first orbital
%       orbsym2   - Symbolic representation of second orbital
%
%   Outputs:
%       TBSK_hop  - Symbolic Slater-Koster hopping term
%       Coff      - Coefficients for sigma and pi bonding terms
%
%   Note:
%       Automatically handles p-s orbital swapping with sign change
%
%   See also HR.TBSK_Coeff_gen, HR.delta_orb, HR.delta_orb_sym
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
