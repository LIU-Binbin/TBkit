function [Rnn_list,nn_smart_t] = nn_smart_t_gen(Atom_smart_t,Rm,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut)
nn_smart_t.seq_from = Atom_smart_t.seq_from;
nn_smart_t.seq_to = Atom_smart_t.seq_to;
count = 1;
reducible_num=(2*search_rangex+1)*(2*search_rangey+1)*(2*search_rangez+1);
Rnn_list = zeros(reducible_num,1);
nn_t = struct('R_vector',[],'R_fractional_diff',[],'Rlength',[],'nn_level',[],'hop_pre',[],'hop',[]);
nn = repmat(nn_t,[reducible_num 1]);
for Rf_a1=-search_rangex:search_rangex
for Rf_a2=-search_rangey:search_rangey
for Rf_a3=-search_rangez:search_rangez
R_vector = [Rf_a1 Rf_a2 Rf_a3];
Rij_cart = (R_vector + Atom_smart_t.R_fractional_diff)*Rm ;
Rlength = norm(Rij_cart);
Rlength = roundn(Rlength,-Accuracy);
Rlmn = Rij_cart/Rlength;
if  0 < Rlength && Rlength < Rlength_cut
orb1 = Atom_smart_t.l_name_from ;
orb2 = Atom_smart_t.l_name_to   ;
orb_sym1 =  Atom_smart_t.orb_sym_from;
orb_sym2 =  Atom_smart_t.orb_sym_to;
orbsym1_n =  HR.subs_xyz(orb_sym1 ,Rlmn);
orbsym2_n =  HR.subs_xyz(orb_sym2 ,Rlmn);
TBSK_hop = HR.TBSK_hop_gen(orb1,orb2,orbsym1_n,orbsym2_n,orb_sym1,orb_sym2);
nn(count,1).hop_pre = TBSK_hop;
nn(count,1).R_vector = R_vector;
nn(count,1).R_fractional_diff = Atom_smart_t.R_fractional_diff;
nn(count,1).Rlength = Rlength;
Rnn_list(count,:) = Rlength;
count = count +1;
end
end
end
end
if count <= reducible_num
Rnn_list(count:reducible_num,:) = [];
nn(count:reducible_num,:) = [];
end
nn_smart_t.nn = nn;
end
