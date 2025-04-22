function H_hr = plus(A,B)
% PLUS Overloaded addition operator for HR objects
%
%   H_HR = PLUS(A,B) implements matrix addition for HR objects
%
%   Inputs:
%       A - HR object or numeric matrix
%       B - HR object or numeric matrix
%   Output:
%       H_hr - Resulting HR object after addition
%
%   Notes:
%       - Handles three cases: HR+HR, HR+numeric, numeric+HR
%       - Checks WAN_NUM compatibility for HR+HR case
%       - Supports both numeric and symbolic operations
%       - Preserves vector hopping information when present
if isa(A,'HR') && isa(B,'HR')
    H_hr1 = A;
    H_hr2 = B;
    if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
        error('WAN_NUM different');
    end
    if ~strcmp(A.Type,'list') && ~strcmp(B.Type,'list')
        H_hr = H_hr1;
        for i = 1:H_hr2.NRPTS
            vector = H_hr2.vectorL(i,:);
            if H_hr.coe
                amp_sym = H_hr2.HcoeL(:,:,i);
                H_hr = H_hr.set_hop_mat(amp_sym,vector,'symadd');
            end
            if H_hr.num
                amp = H_hr2.HnumL(:,:,i);
                H_hr = H_hr.set_hop_mat(amp,vector,'add');
            end
        end
    elseif strcmp(H_hr1.Type,'list') && strcmp(H_hr2.Type,'list')
        [vectorList,~,~] = unique([H_hr1.vectorL;H_hr2.vectorL],'rows');
        C = setdiff(vectorList,H_hr1.vectorL,'rows');
        D = setdiff(vectorList,H_hr2.vectorL,'rows');
        H_hr1 = H_hr1.add_empty_one(C);
        H_hr2 = H_hr2.add_empty_one(D);
        [~,seqA] = ismember(vectorList,H_hr1.vectorL,'rows');
        [~,seqB] = ismember(vectorList,H_hr2.vectorL,'rows');
        H_hr = H_hr1.reseq(':',seqA);
        if H_hr.vectorhopping && B.vectorhopping
            H_hr.AvectorL = H_hr.AvectorL + H_hr2.AvectorL(seqB,:);
            H_hr.BvectorL = H_hr.BvectorL + H_hr2.BvectorL(seqB,:);
            CL1 = H_hr2.CvectorL(1:end/2,:);
            CL2 = H_hr2.CvectorL(end/2+1:end,:);
            H_hr.CvectorL = H_hr.CvectorL + [CL1(seqB,:);CL2(seqB,:)];
            return;
        end
        if H_hr.coe
            H_hr.HcoeL = H_hr.HcoeL + H_hr2.HcoeL(seqB) ;
        end
        if H_hr.num
            H_hr.HnumL = H_hr.HnumL + H_hr2.HnumL(seqB) ;
        end
    elseif ~strcmp(H_hr1.Type,H_hr2.Type)
        H_hr1 = H_hr1.rewrite();
        H_hr2 = H_hr2.rewrite();
        H_hr = H_hr1+H_hr2;
        return;
    end
elseif isa(A,'HR') && ~isa(B,'HR')
    H_hr1 = A;
    H_hr = H_hr1;
    i = H_hr.Line_000;
    if isa(B,'sym')
        H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i) + B;
    elseif isa(B,'Term')
        H_hk_tmp = HK(H_hr.WAN_NUM,2);
        H_hk_tmp.Rm = H_hr.Rm;
        H_hk_tmp = H_hk_tmp + B;
        H_hr = H_hr.plus(H_hk_tmp.kp2TB());
    elseif isa(B,'numeric')
        H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i) + B;
    else
    end
elseif ~isa(A,'HR') && isa(B,'HR')
    H_hr1 = B;
    H_hr = H_hr1;
    if isa(A,'sym')
        for i = 1:H_hr.NRPTS
            H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i) + A;
        end
    elseif isa(A,'numeric')
        for i = 1:H_hr.NRPTS
            H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i) + A;
        end
    else
    end
end
end
