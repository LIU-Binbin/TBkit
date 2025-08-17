clear
sigma_0 = [1 0;0 1];
sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0];
sigma_z = [1 0;0 -1];

Graphene = HR(2);
Graphene = Graphene<'POSCAR';
Graphene = Graphene<'KPOINTS';
search_range = [1 1 0];
maxR = 2.5;
Accuracy = 1e-3;
Graphene = Graphene.nn(search_range, Accuracy ,maxR);
[Rnn,~,~,~] = Graphene.nn_information();
%% 
% 我们指定level_cut = 1，也就是最高只考虑到最近邻项，跨原胞的跃迁也仅考虑到相邻原胞，传入函数，得到的是一个符号化的TB模型，符号化的哈密顿量矩阵储存在HR类的HcoeL属性中。调用symvar_list可以显示自动生成的待定参数名称。此时，其HnumL属性是由空矩阵构成的，手动为生成参数赋值后，就可以用Subsall转化得到数值化的TB模型。查看其下的HnumL，此时应当已成为了非空矩阵，求解本征值后，就可以在布里渊区绘制能带。

Graphene = Graphene.H_TBSK_gen('level_cut',1,'per_dir',[1 1 0]);
Graphene.symvar_list
VppP_1 = 1;
VppP_2 = 0.2;
level_cut = 1;
% list mode
Graphene_list = rewrite(Graphene);
Graphene = rewrite(Graphene_list,'rewind',true);
Graphene_list.symvar_list;
list(Graphene_list);
printout(Graphene);
%Graphene_list.show('HOPPING','scale', 2.4560000896,'atomscale',0.1,'TwoD',true);
%%
Graphene_hk = Graphene_list.HR2HK();
sym(Graphene_hk)
%%
Graphene_htrig = Graphene_list.HR2Htrig();
latex(Graphene_htrig);
%%
K1 = [1/3,1/3,0];
Graphene_hk_K1 = Graphene_htrig.Htrig2HK(K1,'sym',true,'Order',1);
%%
K2 = [-1/3,-1/3,0];
Graphene_hk_K2 = Graphene_htrig.Htrig2HK(K2,'sym',true,'Order',1);
%%
G = [0,0,0];
Graphene_hk_G= Graphene_htrig.Htrig2HK(G,'sym',true,'Order',2);
%%
% C= simplify(sym(Graphene_htrig))
% simplify(C(2,1))
% 
% syms VppP_1 v_k k_x k_y k_plus real
% A= simplify(sym(Graphene_hk_K1))
% A = simplify(subs(A,[VppP_1] ,[2*v_k/(sqrt(3))] ))
% S = sym([exp(1i*pi/3) 0;0 exp(-1i*pi/3)])
% simplify(S*A*S' )
% 
% B = simplify(sym(Graphene_hk_K2))
% B = simplify(subs(B,[VppP_1] ,[2*v_k/(sqrt(3))] ))
% simplify(S'*B*S )
% 
% syms k phi_k real;
% A_prime =  simplify( rewrite(subs(A,[k_x k_y],[k*cos(phi_k) k*sin(phi_k)]),'exp'));
% A_prime = simplify(subs(A_prime,phi_k,phi_k+pi/3))
% 
% B_prime =  simplify( rewrite(subs(B,[k_x k_y],[k*cos(phi_k) k*sin(phi_k)]),'exp'));
% B_prime = simplify(subs(B_prime,phi_k,phi_k+pi/3))
% 
% D = subs(eig(S*A*S'),v_k,1);
% D(1) = -sqrt(D(1)^2);
% D(2) = sqrt(D(2)^2);
% fsurf(D)
% 
% syms VppP_1 v_k k_x k_y k_plus real
% H= simplify(sym(Graphene_hk_K3))
% H = simplify(subs(H,[VppP_1] ,[2*v_k/(sqrt(3))] ))
% S = sym([exp(1i*pi/6) 0;0 exp(-1i*pi/6)])
% simplify(S'*H*S)
