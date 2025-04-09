function H_hk = init(H_hk)
%INIT Initialize Hamiltonian coefficient matrices
%
% Syntax:
%   H_hk = init(H_hk)
%
% Input:
%   H_hk - HK object to initialize
%
% Output:
%   H_hk - HK object with initialized coefficient matrices
%
% Description:
%   Initializes the symbolic coefficient matrices (HcoeL) with:
%   - Real parts 'A_ij' and imaginary parts 'B_ij'
%   - Enforced Hermitian structure (upper triangular + conjugate transpose)
%
% Note:
%   Used when creating new Hamiltonians before term assignment
%   Preserves existing matrix dimensions and Kinds count
%
% Example:
%   Hk = Hk.init(); % Prepare Hamiltonian for term assignment)
sizeH = size(H_hk.HcoeL);
HcoeL_tmp = sym('A',sizeH,'real')+1i*sym('B',sizeH,'real');
for i =1:H_hk.Kinds
    HcoeL_tmp(:,:,i) = triu(HcoeL_tmp(:,:,i))+triu(HcoeL_tmp(:,:,i),1)';
end
H_hk.HcoeL = HcoeL_tmp;
end
