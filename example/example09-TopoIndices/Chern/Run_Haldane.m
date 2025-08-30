sigma_0 = pauli_matrix(0); sigma_x = pauli_matrix(1); sigma_y = pauli_matrix(2); sigma_z = pauli_matrix(3);

Graphene = HR(2);
Graphene = Graphene < 'POSCAR';
Graphene = Graphene < 'KPOINTS';
search_range = [1 1 0];
maxR = 2.5;
Accuracy = 1e-6;
Graphene = Graphene.nn(search_range, Accuracy ,maxR);
[Rnn,~,~,~] = Graphene.nn_information();
%% 
Graphene = Graphene.H_TBSK_gen('level_cut',1,'per_dir',[1 1 0]);
VppP_1 = 1;
%%
V_NNN = [1 0 0; 0 1 0; -1 -1 0];
VppP_2 = 0.1 * VppP_1;

Haldane = Graphene;

Haldane = Haldane.set_hop( -1i*VppP_2,1,1,V_NNN,'sym');
Haldane = Haldane.set_hop(  1i*VppP_2,2,2,V_NNN,'sym');
 
Haldane = Haldane.autohermi();
Haldane_n = Haldane.Subsall();
%%
Chern = Haldane_n.ChernNumber('BAND_index',1);
