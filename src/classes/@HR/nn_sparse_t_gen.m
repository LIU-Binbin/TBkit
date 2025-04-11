function [nn_sparse_temp,Rnn_list] = nn_sparse_t_gen(site1,site2,Rm,search_rangex,search_rangey,search_rangez,Accuracy,Rlength_cut)
% NN_SPARSE_T_GEN Generate sparse neighbor table entry
%
%   [NN_SPARSE_TEMP,RNN_LIST] = NN_SPARSE_T_GEN(SITE1,SITE2,RM,...)
%   generates a single sparse neighbor table entry
%
%   Inputs:
%       site1, site2 - Site information structures
%       Rm - Lattice vectors
%       search_rangex/y/z - Search ranges in each dimension
%       Accuracy - Rounding accuracy
%       Rlength_cut - Cutoff distance
%   Outputs:
%       nn_sparse_temp - Sparse matrix entry data
%       Rnn_list - List of neighbor distances
%
%   Notes:
%       - Helper function for nn_sk_sparse
%       - Uses compact numeric storage format
%       - Filters neighbors by cutoff distance
Rc1 = [site1.rc1,site1.rc2,site1.rc3];
Rc2 = [site2.rc1,site2.rc2,site2.rc3];
R_fractional_diff = -(Rc1 - Rc2);
count = 1;
reducible_num=(2*search_rangex+1)*(2*search_rangey+1)*(2*search_rangez+1);
Rnn_list = zeros(reducible_num,1);
nn_sparse_temp= zeros(reducible_num,10);
for Rf_a1=-search_rangex:search_rangex
for Rf_a2=-search_rangey:search_rangey
for Rf_a3=-search_rangez:search_rangez
R_vector = [Rf_a1 Rf_a2 Rf_a3];
Rij_cart = (R_vector + R_fractional_diff)*Rm ;
Rlength = norm(Rij_cart);
Rlength = roundn(Rlength,-Accuracy);
Rlmn = Rij_cart/Rlength;
if  0 < Rlength && Rlength < Rlength_cut
orb1 = site1.orb   ;
orb2 = site2.orb   ;
orb_sym1 =  site1.orb_sym;
orb_sym2 =  site2.orb_sym;
orbsym1_n =  HR.subs_xyz(orb_sym1 ,Rlmn);
orbsym2_n =  HR.subs_xyz(orb_sym2 ,Rlmn);
nn_sparse_temp(count,3) = Rf_a1;
nn_sparse_temp(count,4) = Rf_a2;
nn_sparse_temp(count,5) = Rf_a3;
nn_sparse_temp(count,6) =  Rlength;
[~,Coff] = HR.TBSK_hop_gen(orb1,orb2,orbsym1_n,orbsym2_n,orb_sym1,orb_sym2);
nn_sparse_temp(count,8) =   Coff(1);
nn_sparse_temp(count,9) =   Coff(2);
nn_sparse_temp(count,10) =  Coff(3);
Rnn_list(count,:) = Rlength;
count = count +1;
end
end
end
end
if count <= reducible_num
Rnn_list(count:reducible_num,:) = [];
nn_sparse_temp(count:reducible_num,:) = [];
end
end
