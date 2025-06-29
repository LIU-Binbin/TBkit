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
%% 必须指明晶系
Msg = read_Magnetic_Sym_El([152,34],"hexagonal");
gen_list = Msg;
%% 3rd 
jahn_symbol_Str = 'aeVVV';
Tensor = jahn_symbol(jahn_symbol_Str);

for i = 1:length(gen_list)
    Tensor = group_transformation(Tensor, gen_list(i));
end
% name of the numerical tensor
[~, ~, SymMatDisplay] = pretty(Tensor,"Table");
SymMatDisplay
%%
jahn_symbol_Str = 'eV4';
Tensor4 = jahn_symbol(jahn_symbol_Str);

for i = 1:length(gen_list)
    Tensor4 = group_transformation(Tensor4, gen_list(i));
end
% name of the numerical tensor
[Table, SymMat, SymMatDisplay]= pretty(Tensor4,"Table");
SymMatDisplay