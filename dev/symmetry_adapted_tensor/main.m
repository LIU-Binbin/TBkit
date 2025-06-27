%% Symmetry Adapted Tensor
%% 手动定义对称性操作的例子

% P = Oper.inversion();
% Mx = Oper.mirror([1,0,0]);
% My = Oper.mirror([0,1,0]);
% Mz = Oper.mirror([0,0,1]);
% C6x = Oper.rotation(1/6,[1,0,0]);
% C6z = Oper.rotation(1/6,[0,0,1]);
% C3x = Oper.rotation(1/3,[1,0,0]);
% C3z = Oper.rotation(1/3,[0,0,1]);
% S6x = Mx*C6x;
% S6z = Mz*C6z;
%% 从Mvasp2trace的磁群对称性表中统一读入对称性操作，必须指明晶系

%msg_BNS_number = "194_268";
%Msg = read_Magnetic_Sym_El("Magnetic_Sym_El/Magnetic_Sym_El_"+msg_BNS_number+".txt","hexagonal");

Msg = read_Magnetic_Sym_El([152,34],"cubic");

gen_list = Msg([1,3,4,7])
% %% basic setting of the tensor,
% tensor_rank = 3;
% is_pesudo_tensor = true;
% is_Time_reversal_odd = true;
% %% initialize the unsymmed tensor
% n_ele = 3^tensor_rank;
% tensor0 = eye(n_ele);
% tensor0 = reshape(tensor0, [ones(1,tensor_rank)*3, n_ele]);
%% 3rd 
tic;
jahn_symbol_Str = 'aeVVV';
Tensor = jahn_symbol(jahn_symbol_Str);
%%
for i = 1:length(gen_list)
    Tensor = group_transformation(Tensor, gen_list(i));
end
%% name of the numerical tensor

[Table,SymMat] = pretty(Tensor,"Table")
toc;
%% 对于五阶自旋霍尔，筛选面内驱动的面外极化面外流

tic;
jahn_symbol_Str = 'eV4';
Tensor5 = jahn_symbol(jahn_symbol_Str);
%%
for i = 1:length(gen_list)
    Tensor5 = group_transformation(Tensor5, gen_list(i));
end
% name of the numerical tensor
[Table,SymMat,SymMatDisplay]= pretty(Tensor5,"Table");
SymMatDisplay
toc;
% bug here!;
% https://www.cryst.ehu.es/cgi-bin/cryst/programs/mtensor_build.pl?magnum=152.34&generators=&transfmat=&jahnsymbol=eV4&database_type=mtensor&trm=a,b,c&sym_cards=&input_type=&x1=1&x2=0&x3=0&x4=0&y1=0&y2=1&y3=0&y4=0&z1=0&z2=0&z3=1&z4=0&pointo=
return;
%% Test for {V}
% bug here ! 
% check difference
% https://www.cryst.ehu.es/cgi-bin/cryst/programs/mtensor_build.pl?magnum=194.268&generators=&transfmat=&jahnsymbol=a{V2}VV&database_type=mtensor&trm=a,b,c&sym_cards=&input_type=&x1=1&x2=0&x3=0&x4=0&y1=0&y2=1&y3=0&y4=0&z1=0&z2=0&z3=1&z4=0&pointo=
tic;
jahn_symbol_Str = 'a{V2}VV';
Tensor4 = jahn_symbol(jahn_symbol_Str);
%
for i = 1:length(gen_list)
    Tensor4 = group_transformation(Tensor4, gen_list(i));
end
% name of the numerical tensor

[Table,SymMat,SymMatDisplay]= pretty(Tensor4,"Table");
SymMatDisplay
toc;