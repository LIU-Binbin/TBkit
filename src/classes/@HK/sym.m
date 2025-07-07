function H_hk_sym = sym(H_hk)
%SYM Convert to symbolic matrix
%
% Syntax:
%   H_sym = sym(Hk)
%
% Input:
%   Hk - HK Hamiltonian
%
% Output:
%   H_sym - Symbolic matrix representation
%
% Description:
%   Extracts full symbolic Hamiltonian (Hk_sym property)
%   Equivalent to Hk.Hk_sym but with standard conversion interface
%
% Example:
%   H_mat = sym(Hk); % Get symbolic matrix
if H_hk.num
    H_hk_sym = sym(zeros(H_hk.Basis_num,H_hk.Basis_num));
    for i =1:H_hk.Kinds
        H_hk_sym = H_hk_sym + H_hk.HnumL(:,:,i)*H_hk.HsymL_k(i);
    end
else
    H_hk_sym = H_hk.Hk_sym;
end
end
