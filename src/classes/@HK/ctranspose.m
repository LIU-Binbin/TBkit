function H_hk = ctranspose(H_hk)
%CTRANSPOSE Hermitian conjugate of HK Hamiltonian
%
% Syntax:
%   H_hk = ctranspose(H_hk)
%   H_hk = H_hk'
%
% Input:
%   H_hk - HK object to transpose
%
% Output:
%   H_hk - Hermitian conjugate of Hamiltonian
%
% Description:
%   Computes the Hermitian conjugate (conjugate transpose) of the Hamiltonian
%   by transposing both numeric (HnumL) and symbolic (HcoeL) coefficient
%   matrices while applying complex conjugation.
%
% Note:
%   Implements the ' operator for HK class objects
for i =1:H_hk.Kinds
    H_hk.HnumL(:,:,i) = H_hk.HnumL(:,:,i)';
    H_hk.HcoeL(:,:,i) = H_hk.HcoeL(:,:,i)';
end
end
