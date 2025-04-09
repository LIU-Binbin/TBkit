function C = horzcat(A,B)
%HORZCAT Horizontal concatenation of HK Hamiltonians
%
% Syntax:
%   C = [A B]       % Operator form
%   C = horzcat(A,B) % Functional form
%
% Inputs:
%   A, B - HK objects to concatenate
%
% Output:
%   C - Combined HK Hamiltonian with block-diagonal structure
%
% Description:
%   Performs horizontal concatenation of two HK objects by:
%   1. Combining their basis sets
%   2. Creating block-diagonal matrices for:
%      - Numeric coefficients (HnumL)
%      - Symbolic coefficients (HcoeL)
%   3. Preserving all term structure (Kinds, Degree)
%
% Note:
%   Only implemented for HK-HK concatenation
%   Resets Trig_to_save property to zeros
%
% Example:
%   H_combined = [Hk1 Hk2]; % Creates block-diagonal Hamiltonian
if isa(A,'HK') && isa(B,'HK')
    H_hk1 = A;
    H_hk2 = B;
    H_hk = A;
    H_hk.Basis = [H_hk1.Basis;H_hk2.Basis];
    H_hk.Basis_num = H_hk1.Basis_num+H_hk2.Basis_num;
    H_hk.HnumL = zeros(H_hk.Basis_num,H_hk.Basis_num,H_hk.Kinds);
    H_hk.HcoeL = sym(H_hk.HnumL);
    for i = 1:H_hk.Kinds
        H_hk.HnumL(:,:,i) = blkdiag(H_hk1.HnumL(:,:,i) ,...
            H_hk2.HnumL(:,:,i));
        H_hk.HcoeL(:,:,i) = blkdiag(H_hk1.HcoeL(:,:,i) ,...
            H_hk2.HcoeL(:,:,i));
    end
    H_hk.Trig_to_save =sym(zeros(H_hk.Basis_num,H_hk.Basis_num));
else
end
C = H_hk;
end
