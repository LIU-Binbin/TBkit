%% HRclass Tutorial
% 基于HRclass实现石墨烯model      
% 
% 2021.08
%% 
% * Author: parkman
% * Email：parkman@buaa.edu.cn
%% 

syms t real;
Grpahene_TB = HR(2);
A_vectorL = [0,0,0;1,0,0;0,-1,0];
B_vectorL = [0,0,0;-1,0,0;0,1,0];
%%
Grpahene_TB = Grpahene_TB.set_hop(t,1,2,A_vectorL,'sym');
Grpahene_TB = Grpahene_TB.set_hop(t,2,1,B_vectorL,'sym')
%%
Grpahene_TB.printout;
%%
Grpahene_TB_list = Grpahene_TB.rewrite()

%%
Grpahene_TB_list.list();
%%
Rm = [1,0,0;-0.5,sqrt(3)/2,0;0,0,1];
Grpahene_TB_list = Grpahene_TB_list.input_Rm(Rm);
orbL = [2/3,1/3,0;1/3,2/3,0];
Grpahene_TB_list.orbL = orbL;
%%
Grpahene_TB_list.show('HOPPING','scale', 2.4560000896,'atomscale',1,'TwoD',true);
%% 补充用纯代码 作k-path的后续内容