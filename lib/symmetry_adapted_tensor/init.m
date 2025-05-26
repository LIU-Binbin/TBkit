function Tensor = init(Tensor)
tensor_rank = Tensor.tensor_rank;
n_ele = 3^tensor_rank;
tensor_symmed_mat = eye(n_ele);
%
tensor_names = sym('chi_', ones(1,tensor_rank)*3, 'real');
tensor_names_vector = reshape(tensor_names, n_ele, 1);
tensor_names_vector_Str = string(tensor_names_vector);
tensor_names_vector_StrL  = split(tensor_names_vector_Str,{'chi_','_'});
tensor_names_vector_StrL(:,1) =[];
Num_Mat = double(tensor_names_vector_StrL);

%
Npair = size(Tensor.SymmetricPair,1);
for i = 1:Npair
    SymmetricPair = Tensor.SymmetricPair(i,:);
    %%
    Ni = SymmetricPair(1);
    Nj = SymmetricPair(2);

    Nselect = (Num_Mat(:,Ni) == Ni & Num_Mat(:,Nj) == Nj) | ...
        (Num_Mat(:,Nj) == Ni & Num_Mat(:,Ni) == Nj) ;
    Nselect = 1:n_ele;
    OriginSeq = find(Nselect);
    Num_Mat_Sel = Num_Mat(Nselect,:);
    Num_Mat_Ex = Num_Mat_Sel;
    Num_Mat_Ex(:,Ni) = Num_Mat_Sel(:,Nj);
    Num_Mat_Ex(:,Nj) = Num_Mat_Sel(:,Ni);

    [~,CheckList] = ismember(Num_Mat_Ex,Num_Mat_Sel,'rows');
    % tensor_symmed_mat(Ni,:) = tensor_symmed_mat(Ni,:) + tensor_symmed_mat(Nj,:) ;
    % tensor_symmed_mat(Nj,:) = tensor_symmed_mat(Nj,:) + tensor_symmed_mat(Ni,:) ;
    for j = 1:size(CheckList,1)
        Origin_Ni = OriginSeq(j);
        Origin_Nj = OriginSeq(CheckList(j));
        tensor_symmed_mat(Origin_Ni,:) = tensor_symmed_mat(Origin_Ni,:);
        tensor_symmed_mat(Origin_Nj,:) = tensor_symmed_mat(Origin_Ni,:);
    end
end
%
Npair = size(Tensor.AsymmetricPair,1);
for i = 1:Npair
    SymmetricPair = Tensor.AsymmetricPair(i,:);
    %%
    Ni = SymmetricPair(1);
    Nj = SymmetricPair(2);

    Nselect = (Num_Mat(:,Ni) == Ni & Num_Mat(:,Nj) == Nj) | ...
        (Num_Mat(:,Nj) == Ni & Num_Mat(:,Ni) == Nj) ;
    Nselect = 1:n_ele;
    OriginSeq = find(Nselect);
    Num_Mat_Sel = Num_Mat(Nselect,:);
    Num_Mat_Ex = Num_Mat_Sel;
    Num_Mat_Ex(:,Ni) = Num_Mat_Sel(:,Nj);
    Num_Mat_Ex(:,Nj) = Num_Mat_Sel(:,Ni);

    [~,CheckList] = ismember(Num_Mat_Ex,Num_Mat_Sel,'rows');
    % tensor_symmed_mat(Ni,:) = tensor_symmed_mat(Ni,:) + tensor_symmed_mat(Nj,:) ;
    % tensor_symmed_mat(Nj,:) = tensor_symmed_mat(Nj,:) + tensor_symmed_mat(Ni,:) ;
    for j = 1:size(CheckList,1)
        Origin_Ni = OriginSeq(j);
        Origin_Nj = OriginSeq(CheckList(j));
        if Origin_Ni ~= Origin_Nj
            tensor_symmed_mat(Origin_Ni,:) = tensor_symmed_mat(Origin_Ni,:)-tensor_symmed_mat(Origin_Nj,:);
            tensor_symmed_mat(Origin_Nj,:) = tensor_symmed_mat(Origin_Nj,:)-tensor_symmed_mat(Origin_Ni,:);
        else
            tensor_symmed_mat(Origin_Ni,:) = tensor_symmed_mat(Origin_Ni,:)-tensor_symmed_mat(Origin_Ni,:);
            tensor_symmed_mat(Origin_Nj,:) = tensor_symmed_mat(Origin_Nj,:)-tensor_symmed_mat(Origin_Nj,:);
        end
    end
end
%
Npair = size(Tensor.TsymmetricPair,1);
for i = 1:Npair
    SymmetricPair = Tensor.TsymmetricPair(i,:);
    %%
    Ni = SymmetricPair(1);
    Nj = SymmetricPair(2);

    Nselect = (Num_Mat(:,Ni) == Ni & Num_Mat(:,Nj) == Nj) | ...
        (Num_Mat(:,Nj) == Ni & Num_Mat(:,Ni) == Nj) ;
    Nselect = 1:n_ele;
    OriginSeq = find(Nselect);
    Num_Mat_Sel = Num_Mat(Nselect,:);
    Num_Mat_Ex = Num_Mat_Sel;
    Num_Mat_Ex(:,Ni) = Num_Mat_Sel(:,Nj);
    Num_Mat_Ex(:,Nj) = Num_Mat_Sel(:,Ni);

    [~,CheckList] = ismember(Num_Mat_Ex,Num_Mat_Sel,'rows');
    % tensor_symmed_mat(Ni,:) = tensor_symmed_mat(Ni,:) + tensor_symmed_mat(Nj,:) ;
    % tensor_symmed_mat(Nj,:) = tensor_symmed_mat(Nj,:) + tensor_symmed_mat(Ni,:) ;
    for j = 1:size(CheckList,1)
        Origin_Ni = OriginSeq(j);
        Origin_Nj = OriginSeq(CheckList(j));
        tensor_symmed_mat(Origin_Ni,:) = tensor_symmed_mat(Origin_Ni,:);
        tensor_symmed_mat(Origin_Nj,:) = tensor_symmed_mat(Origin_Ni,:);
    end
end
%
Npair = size(Tensor.TasymmetricPair,1);
for i = 1:Npair
    SymmetricPair = Tensor.TasymmetricPair(i,:);
    %%
    Ni = SymmetricPair(1);
    Nj = SymmetricPair(2);

    Nselect = (Num_Mat(:,Ni) == Ni & Num_Mat(:,Nj) == Nj) | ...
        (Num_Mat(:,Nj) == Ni & Num_Mat(:,Ni) == Nj) ;
    Nselect = 1:n_ele;
    OriginSeq = find(Nselect);
    Num_Mat_Sel = Num_Mat(Nselect,:);
    Num_Mat_Ex = Num_Mat_Sel;
    Num_Mat_Ex(:,Ni) = Num_Mat_Sel(:,Nj);
    Num_Mat_Ex(:,Nj) = Num_Mat_Sel(:,Ni);

    [~,CheckList] = ismember(Num_Mat_Ex,Num_Mat_Sel,'rows');
    % tensor_symmed_mat(Ni,:) = tensor_symmed_mat(Ni,:) + tensor_symmed_mat(Nj,:) ;
    % tensor_symmed_mat(Nj,:) = tensor_symmed_mat(Nj,:) + tensor_symmed_mat(Ni,:) ;
    for j = 1:size(CheckList,1)
        Origin_Ni = OriginSeq(j);
        Origin_Nj = OriginSeq(CheckList(j));
        tensor_symmed_mat(Origin_Ni,:) = tensor_symmed_mat(Origin_Ni,:);
        tensor_symmed_mat(Origin_Nj,:) = -tensor_symmed_mat(Origin_Ni,:);
    end
end

Tensor.tensor0 = reshape(tensor_symmed_mat, [ones(1,tensor_rank)*3, n_ele]);
% tensor0 =
end