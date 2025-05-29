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
msg_BNS_number = "194_268";
gen_list = read_Magnetic_Sym_El("Magnetic_Sym_El/Magnetic_Sym_El_"+msg_BNS_number+".txt","hexagonal"); 
%% basic setting of the tensor
tensor_rank = 5;
is_pesudo_tensor = true;
is_Time_reversal_odd = true;
%% initialize the unsymmed tensor
n_ele = 3^tensor_rank;
tensor0 = eye(n_ele);
tensor0 = reshape(tensor0, [ones(1,tensor_rank)*3, n_ele]);
%%
for i = 1:length(gen_list)
    tensor_symmed = group_transformation(tensor0, gen_list(i), ...
        "is_pesudo_tensor",is_pesudo_tensor, ...
        "is_Time_reversal_odd",is_Time_reversal_odd);

    tensor_symmed_mat = reshape(tensor_symmed, n_ele, n_ele);
    [tensor_symmed_reduce, independent_ele] = rref(tensor_symmed_mat');

    tensor0 = reshape(tensor_symmed_reduce', [ones(1,tensor_rank)*3, n_ele]);
end

%% name of the numerical tensor
tensor_names = sym('X_', ones(1,tensor_rank)*3, 'real');

% 对于五阶自旋霍尔，筛选面内驱动的面外极化面外流
% tensor_names_tmp = tensor_names;
% tensor_names_tmp(3,3,1:2,1:2,1:2) = 0;
% tensor_names = tensor_names - tensor_names_tmp;

tensor_names_vector = reshape(tensor_names, n_ele, 1);
n_indp = length(independent_ele);

for i = 1:n_ele
    % 增加一步筛选
    if ~all(tensor_symmed_reduce(1:n_indp,i)==0) && tensor_names_vector(i) ~=0
        eq = sum(tensor_symmed_reduce(1:n_indp,i) .* tensor_names_vector(independent_ele));
        disp(string(tensor_names_vector(i))+ ...
            " = "+ ...
            string(eq))
    end
end