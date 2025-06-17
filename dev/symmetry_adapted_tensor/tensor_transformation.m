function tensor1 = tensor_transformation(Tensor, OperObj)
arguments
    Tensor 
    OperObj Oper
end
% n_ele = Tensor.n_ele;
tensor_rank = Tensor.tensor_rank;
tensor0 = Tensor.tensor0;


OperR = double(OperObj.R);

tensor1 = tensor0;
for i = 1:tensor_rank
    % dh*efghL -> defgL
    tensor1 = tensorprod(OperR, tensor1, 2, tensor_rank);

    % if class(tensor0) == "sym"
    %     tensor1 = reshape(tensor1, 3, []);
    % 
    %     % ae*efghL -> afghL, abcdefgh are the basis, L is unchanged index
    %     tensor1 = operR * tensor1;
    %     tensor1 = reshape(tensor1, [3*ones(1,tensor_rank),3^tensor_rank]);
    % 
    %     % afghL -> fghaL
    %     tensor1 = permute(tensor1, [2:tensor_rank, 1, tensor_rank+1]);
    % elseif class(tensor0) == "double"
    %     operR = double(operR);
    %     % dh*efghL -> defgL
    %     tensor1 = tensorprod(operR, tensor1, 2, tensor_rank);
    % end  
end

if Tensor.is_pesudo_tensor
    tensor1 = tensor1 * det(OperR);
end

if Tensor.is_Time_reversal_odd && OperObj.conjugate
    tensor1 = tensor1 * -1;
end

end