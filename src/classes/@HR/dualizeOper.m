function [H_hr, VectorDistMat] = dualizeOper(H_hr, SymOper,options)
arguments
    H_hr 
    SymOper 
    options.associate_orbL = [];
    options.eta = 1e-6;
end
% DUALIZEOPER Apply dual operation to Hamiltonian in HR format [Optimized]
eta = options.eta;

Rf = double(Oper.Rc2Rf(inv(SymOper.R), H_hr.Rm));
U = SymOper.U;
NBands = size(U,1);
invRf = inv(Rf);
invU = inv(U);

NonZeroUIndex = find(U~=0);
[Possible_iL,Possible_lL] = ind2sub([NBands,NBands],NonZeroUIndex);


NonZeroInvUIndex = find(invU~=0);
[Possible_mL,Possible_jL] = ind2sub([NBands,NBands], NonZeroInvUIndex);

[valid_T, valid_ilmjL,valid_UinvU] = ...
    Tij2lm(H_hr, Rf,Possible_iL,Possible_lL,Possible_mL,Possible_jL,U,invU,NBands,eta,options.associate_orbL);

% [T,iT] = unique(valid_T,"rows");
% [unique_lmL] = unique(valid_ilmjL(:,2:3),"rows");
% [unique_ijL] = unique(valid_ilmjL(:,[1,4]),"rows");
% 
% 
Dim = H_hr.Dim;
% 

% vectorListUnique = unique(H_hr.vectorL(:, 1:Dim), 'rows');
% 
% m = size(vectorListUnique, 1);
% n = size(iT, 1);
% 
% 
% row_idx = repelem((1:m)', n);       % 行索引扩展
% col_idx = repmat((1:n)', m, 1);    % 列索引扩展
% diffVec = vectorListUnique(row_idx, :) - T(col_idx, :);
% 
% 
% [uniqueRvector,iRvector] = unique(diffVec / Rf, 'rows');
% 
% 
% % diffset(,);
% 
% Add = [gridhorzcat(uniqueRvector,unique_lmL);
%     gridhorzcat(vectorListUnique,unique_ijL)];
% Add = [Add;...
%     -Add(:,1:3),Add(:,[5,4])];
% 
% Add = unique(Add,"rows");
% Add = setdiff(Add,H_hr.vectorL,'rows');
% if ~isempty(Add)
%     H_hr  = H_hr.add_empty_one(Add);
% end

if length(H_hr.vectorL_map) ~= H_hr.NRPTS
    error;
end
NRPTS_ = H_hr.NRPTS;
% vectorList = H_hr.vectorL;
%%
% [sort_valid_ilmjL,sortseq] = sort(valid_ilmjL,[1,4]);
% sort_valid_UinvU =  valid_UinvU(sortseq);
% sort_valid_T = sort_valid_T(sortseq);


% Preallocate based on mode
if H_hr.coe
    VectorDistMat = sparse(NRPTS_, NRPTS_);
    num_mode = 0;
elseif H_hr.num
    VectorDistMat.offset = zeros(NRPTS_, 1,'int32');
    VectorDistMat.idL = [];
    VectorDistMat.UL = [];
    num_mode = 1;
end
intPos = 1;  % 整数数据当前位置
% cmpPos = 1;  % 复数数据当前位置
% ml_pairs  = [];
% Main processing loop
in = 0;
while  in < H_hr.NRPTS
    in = in+1;
    i =  H_hr.vectorL(in,Dim+1);
    j = H_hr.vectorL(in,Dim+2);
    Rvector = H_hr.vectorL(in,1:Dim);
    % Find ml_list for current ij
    Selectseq = valid_ilmjL(:,1) == i & valid_ilmjL(:,4) == j;
    ml_list = valid_ilmjL(Selectseq,2:3);
    U = valid_UinvU(Selectseq);
    Rselect = round((Rvector - valid_T(Selectseq,:)) / Rf);
    %
    CheckList =  cellfun(@mat2str, num2cell([Rselect,ml_list], 2), 'UniformOutput', false) ;
    idj = [];
    for iCheck = 1:numel(CheckList)
        idj_tmp = find(strcmp(CheckList{iCheck},H_hr.vectorL_map));
        if isempty(idj_tmp)
            H_hr = H_hr.add_empty_one(str2num(CheckList{iCheck}));
            idj_tmp = H_hr.NRPTS;
        end
        idj = [idj,idj_tmp];
    end
   
    Nidj = length(idj);
    if length(idj) ~= length(U)
        error
    end

 

    
    
    % Update VectorDistMat based on mode
    if num_mode
        VectorDistMat.offset(in) = intPos;
        if ~isempty(idj)
            VectorDistMat.idL(intPos:(intPos+Nidj-1),1) = idj;
            VectorDistMat.UL(intPos:(intPos+Nidj-1),1) = U;
            intPos = intPos+length(idj);
        else

        end
    else
        % For coe mode: store all connections
        VectorDistMat(in, idj) = U;
    end
    
end

% Nested function for cleaning complex values
    function cleaned = cleanComplex(val, tol)
        real_part = real(val);
        imag_part = imag(val);
        real_part(abs(real_part) < tol) = 0;
        imag_part(abs(imag_part) < tol) = 0;
        cleaned = real_part + 1i * imag_part;
    end
end



function [valid_T, valid_ilmjL,valid_UinvU] = Tij2lm(H_hr, Rf,Possible_iL,Possible_lL,Possible_mL,Possible_jL,U,invU,Nbands,eta,associate_orbL)
arguments
    H_hr 
    Rf 
    Possible_iL 
    Possible_lL 
    Possible_mL 
    Possible_jL 
    U 
    invU 
    Nbands 
    eta = 1e-6;
    associate_orbL = [];
end
% Optimized TIJ2LM with vectorization
% if isempty(associate_orbL)
    orblist = H_hr.orbL;
%     orblist2 = H_hr.orbL;
% else
%     orblist1= H_hr.orbL;
%     orblist2 = H_hr.associate_orbL;
% end
Norb = size(orblist, 1);

Possible_ilL = [Possible_iL,Possible_lL];
Possible_mjL = [Possible_mL,Possible_jL];

Possible_ilmjL = gridhorzcat(Possible_ilL,Possible_mjL);
             
% i l m j
% 1 2 3 4
% 
%T_ij_ml = (m - l)-(j -i) ; 3 2   + 1
T_ij_ml = double( (orblist(Possible_ilmjL(:,3), :) - orblist(Possible_ilmjL(:,2), :) )* Rf + ...
           orblist(Possible_ilmjL(:,1), :) - orblist(Possible_ilmjL(:,4), :) );
lattice_mask = all(abs(T_ij_ml - round(T_ij_ml)) < eta, 2);
valid_ilmjL = Possible_ilmjL(lattice_mask, :);
valid_T     = T_ij_ml(lattice_mask, :);

valid_U  = U(sub2ind([Nbands,Nbands],valid_ilmjL(:,1),valid_ilmjL(:,2)));
valid_invU  = invU(sub2ind([Nbands,Nbands],valid_ilmjL(:,3),valid_ilmjL(:,4)));

valid_UinvU = valid_U.*valid_invU;

end