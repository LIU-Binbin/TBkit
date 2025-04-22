function H_hr = kron(A,B)
%KRON Kronecker tensor product for HR objects
%   Performs Kronecker product between HR objects or between HR and matrices.
%
%   Syntax:
%       H_hr = kron(A,B)
%
%   Inputs:
%       A - First operand (HR object or matrix)
%       B - Second operand (HR object or matrix)
%
%   Outputs:
%       H_hr - Resulting HR object from Kronecker product
%
%   Example:
%       H_kron = kron(H1, H2); % Kronecker product of two HR objects
%       H_kron = kron(H1, eye(2)); % Kronecker product with identity matrix
if isa(A,'HR') && isa(B,'HR')
    H_hr1 = A;
    H_hr2 = B;
    H_hr =  HR(H_hr1.WAN_NUM * H_hr2.WAN_NUM,...
        unique([H_hr1.vectorL;H_hr2.vectorL],'rows'));
    for i = 1:H_hr.NRPTS
        vector = H_hr.vectorL(i,:);
        [~,seq1]=ismember(vector,H_hr1.vectorL,'rows');
        [~,seq2]=ismember(vector,H_hr2.vectorL,'rows');
        if seq1 ~= 0 && seq2 ~=0
            amp = kron(H_hr1.HnumL(:,:,seq1) ,...
                H_hr2.HnumL(:,:,seq2));
            amp_sym = kron(H_hr1.HcoeL(:,:,seq1) ,...
                H_hr2.HcoeL(:,:,seq2));
            H_hr = H_hr.set_hop_mat(amp,vector,'set');
            H_hr = H_hr.set_hop_mat(amp_sym,vector,'sym');
        end
    end
elseif isa(A,'HR') && ~isa(B,'HR')
    H_hr1 = A;
    H_hr = H_hr1;
    N_B = length(B);
    H_hr.orbL = kron(H_hr1.orbL,ones(N_B,1));
    H_hr.quantumL = kron(H_hr1.quantumL,ones(N_B,1));
    H_hr.elementL = kron(H_hr1.elementL,ones(N_B,1));
    H_hr.orb_symL = kron(H_hr1.orb_symL,ones(N_B,1));
    if A.coe
        for i = 1:H_hr.NRPTS
            tmpcoeL(:,:,i) = kron(H_hr.HcoeL(:,:,i),B);
        end
        H_hr.HcoeL = tmpcoeL;
    elseif A.num
        for i = 1:H_hr.NRPTS
            tmpnumL(:,:,i) = kron(H_hr.HnumL(:,:,i),B);
        end
        H_hr.HnumL = tmpnumL;
    else
    end
elseif ~isa(A,'HR') && isa(B,'HR')
    H_hr1 = B;
    H_hr = H_hr1;
    N_A = length(A);
    H_hr.orbL = kron(ones(N_A,1),H_hr1.orbL);
    H_hr.quantumL = kron(ones(N_A,1),H_hr1.quantumL);
    H_hr.elementL = kron(ones(N_A,1),H_hr1.elementL);
    H_hr.orb_symL = kron(ones(N_A,1),H_hr1.orb_symL);
    if B.coe
        for i = 1:H_hr.NRPTS
            tmpcoeL(:,:,i) = kron( A , H_hr.HcoeL(:,:,i));
        end
        H_hr.HcoeL = tmpcoeL;
    elseif B.num
        for i = 1:H_hr.NRPTS
            tmpnumL(:,:,i) = kron( A , H_hr.HnumL(:,:,i));
        end
        H_hr.HnumL = tmpnumL;
    else
    end
end
end
