function TBkitObj = SliceGen(TBkitObj)
%SLICEGEN Generate slicing information for TBkit object
%   TBKITOBJ = SLICEGEN(TBKITOBJ) processes a TBkit object to generate
%   slicing information for sparse matrix operations.
%
%   Input:
%       TBkitObj - Input TBkit object (Htrig or HR type)
%
%   Output:
%       TBkitObj - Modified object with slicing information
%
%   Note:
%       Handles both Htrig and HR object types differently


% $H_{i j}^{\mathbf{k}}=\left\langle\chi_{i}^{\mathbf{k}}|H| \chi_{j}^{\mathbf{k}}\right\rangle=\sum_{\mathbf{R}} e^{i \mathbf{k} \cdot\left(\mathbf{R}+\mathbf{t}_{j}-\mathbf{t}_{i}\right)} H_{i j}(\mathbf{R})$
DIM = TBkitObj.Dim;
if strcmp(TBkitObj.Type,'list')
    switch class(TBkitObj)
        case 'Htrig'
            Kind = size(TBkitObj.HsymL_numL,1);
            [ij_list,index_row] = sortrows(TBkitObj.HsymL_numL(:,DIM+1:DIM+2));
        case 'HR'
            Kind = TBkitObj.NRPTS;
            [ij_list,index_row] = sortrows(TBkitObj.vectorL(:,DIM+1:DIM+2));
    end
    TBkitObj = TBkitObj.reseq(':',index_row);% ?
    [TBkitObj.Sparse_vector,SliceList] = unique(ij_list,'rows');
    TBkitObj.N_Sparse_vector = size(TBkitObj.Sparse_vector,1);
    TBkitObj.CutList = [SliceList,[SliceList(2:end)-1;Kind]];
end
end