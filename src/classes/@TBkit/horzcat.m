function C = horzcat(A,B)
%HORZCAT Horizontal concatenation of HK objects
%   C = horzcat(A, B) combines two HK objects by block diagonal concatenation.
%
%   Inputs:
%       A, B - HK objects to concatenate
%
%   Output:
%       C - Combined HK object with block-diagonal Hamiltonian terms
%
%   Example:
%       H_combined = [H1, H2]; % Uses horzcat internally
if isa(A,'HK') && isa(B,'HK')
    H_hk1 = A;
    H_hk2 = B;
    H_hk = A;
    H_hk.Basis = [H_hk1.Basis; H_hk2.Basis];
    H_hk.Basis_num = H_hk1.Basis_num + H_hk2.Basis_num;
    H_hk.HnumL = zeros(H_hk.Basis_num, H_hk.Basis_num, H_hk.Kinds);
    H_hk.HcoeL = sym(H_hk.HnumL);

    for i = 1:H_hk.Kinds
        H_hk.HnumL(:,:,i) = blkdiag(H_hk1.HnumL(:,:,i), H_hk2.HnumL(:,:,i));
        H_hk.HcoeL(:,:,i) = blkdiag(H_hk1.HcoeL(:,:,i), H_hk2.HcoeL(:,:,i));
    end
    H_hk.Trig_to_save = sym(zeros(H_hk.Basis_num, H_hk.Basis_num));
else
    error('Inputs must be HK objects');
end
C = H_hk;
end
