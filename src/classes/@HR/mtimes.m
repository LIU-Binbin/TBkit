function H_hr = mtimes(A,B)
% MTIMES Overloaded matrix multiplication for HR objects
%
%   H_HR = MTIMES(A,B) implements matrix multiplication (A*B) for HR objects
%
%   Inputs:
%       A - HR object or numeric matrix
%       B - HR object or numeric matrix
%   Output:
%       H_hr - Resulting HR object after multiplication
%
%   Notes:
%       - Handles three cases: HR*HR, HR*numeric, numeric*HR
%       - Checks WAN_NUM compatibility for HR*HR case
%       - Supports both numeric and symbolic operations
if isa(A,'HR') && isa(B,'HR')
    H_hr1 = A;
    H_hr2 = B;
    if H_hr1.WAN_NUM ~= H_hr2.WAN_NUM
        error('WAN_NUM different');
    end
    H_hr =  HR(H_hr2.WAN_NUM,...
        unique([H_hr1.vectorL;H_hr2.vectorL],'rows'));
    for i = 1:H_hr.NRPTS
        vector = H_hr.vectorL(i,:);
        [~,seq1]=ismember(vector,H_hr1.vectorL,'rows');
        [~,seq2]=ismember(vector,H_hr2.vectorL,'rows');
        if seq1 ~= 0 && seq2 ~=0
            amp = H_hr1.HnumL(:,:,seq1) * H_hr2.HnumL(:,:,seq2);
            amp_sym = H_hr1.HcoeL(:,:,seq1) * H_hr2.HcoeL(:,:,seq2);
            H_hr = H_hr.set_hop_mat(amp,vector,'set');
            H_hr = H_hr.set_hop_mat(amp_sym,vector,'sym');
        end
    end
elseif isa(A,'HR') && ~isa(B,'HR')
    H_hr1 = A;
    H_hr = H_hr1;
    if A.coe
        for i = 1:H_hr.NRPTS
            H_hr.HcoeL(:,:,i) = H_hr.HcoeL(:,:,i) * B;
        end
    elseif A.num
        for i = 1:H_hr.NRPTS
            H_hr.HnumL(:,:,i) = H_hr.HnumL(:,:,i) * B;
        end
    else
    end
elseif ~isa(A,'HR') && isa(B,'HR')
    H_hr1 = B;
    H_hr = H_hr1;
    if B.coe
        for i = 1:H_hr.NRPTS
            H_hr.HcoeL(:,:,i) = A * H_hr.HcoeL(:,:,i);
        end
    elseif B.num
        for i = 1:H_hr.NRPTS
            H_hr.HnumL(:,:,i) = A * H_hr.HnumL(:,:,i);
        end
    else
    end
end
end
