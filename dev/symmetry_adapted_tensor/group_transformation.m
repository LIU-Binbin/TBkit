function [Tensor,independent_ele] = group_transformation(Tensor, generator,options)
arguments
    Tensor ;
    generator Oper;
    options.generator = true;
    options.eta = 1e-6;
end
n_ele = Tensor.n_ele;
tensor_rank = Tensor.tensor_rank;
tensor0 = Tensor.tensor0;

if options.generator 
    group = generate_group(generator);
else
    group = generator;
end
Noper = length(group);

tensor_symmed = zeros(size(tensor0));
for j = 1:Noper
    tensor1 = tensor_transformation(Tensor, group(j));
    tensor_symmed = tensor_symmed + tensor1;
end
tensor_symmed = tensor_symmed/Noper;
tensor_symmed( abs(tensor_symmed) < options.eta) = 0;  % remove numerical error
tensor_symmed_mat = reshape(tensor_symmed, n_ele, n_ele);
[tensor_symmed_reduce, independent_ele] = rref(tensor_symmed_mat');

Tensor.tensor0 = reshape(tensor_symmed_reduce', [ones(1,tensor_rank)*3, n_ele]);
end