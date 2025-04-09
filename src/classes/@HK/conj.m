function H_hk = conj(H_hk)
%CONJ Complex conjugate of HK Hamiltonian
%
% Syntax:
%   H_hk = conj(H_hk)
%
% Input:
%   H_hk - HK object to conjugate
%
% Output:
%   H_hk - Hamiltonian with conjugated coefficients
%
% Description:
%   Applies complex conjugation to both numeric (HnumL) and symbolic (HcoeL)
%   coefficient matrices of the Hamiltonian. Preserves the structure while
%   conjugating all elements.
%
% Note:
%   Equivalent to element-wise MATLAB conj() operation for HK class
for i =1:H_hk.Kinds
    H_hk.HnumL(:,:,i) = conj(H_hk.HnumL(:,:,i));
    H_hk.HcoeL(:,:,i) = conj(H_hk.HcoeL(:,:,i));
end
end
