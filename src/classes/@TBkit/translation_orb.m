function [orb_one_incell,translation_vector]=translation_orb(orb_one)
orb_one_incell = mod(orb_one,1);
if orb_one_incell(1) ==1
    orb_one_incell(1)=0;
end
if orb_one_incell(2) ==1
    orb_one_incell(2)=0;
end
if orb_one_incell(3) ==1
    orb_one_incell(3)=0;
end
translation_vector = round(orb_one_incell - orb_one);
end