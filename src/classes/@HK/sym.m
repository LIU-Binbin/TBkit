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
H_hk_sym = H_hk.Hk_sym;
end
