function H_hk = setup_rough(H_hk, symbolic_polynomial, pauli_mat, silence)
%SETUP_ROUGH Configure Hamiltonian from polynomial expression
%
% Syntax:
%   Hk = setup_rough(Hk, poly, mat)
%   Hk = setup_rough(..., silence)
%
% Inputs:
%   Hk - HK Hamiltonian to configure
%   symbolic_polynomial - Polynomial expression
%   pauli_mat - Matrix basis elements
%   silence - Suppress output if true (default=false)
%
% Output:
%   Hk - Configured Hamiltonian
%
% Description:
%   Higher-level setup that:
%   1. Expands and decomposes polynomial
%   2. Extracts coefficients and variables
%   3. Delegates to setup() for actual configuration
%
% Note:
%   Wrapper for more convenient single-term setup
%   Uses standardize_sym for consistent formatting
%
% Example:
%   Hk = Hk.setup_rough(sym('a*k_x'), sigma_x);
if nargin < 4
    silence = false;
end
[coeffs_list,symvar_monomial_list] = coeffs(expand(symbolic_polynomial));
nc = length(coeffs_list);
for i =1:nc
    [Var_cell{i},k_cell{i},~] = HK.coeff_extract(HK.standardize_sym(symvar_monomial_list(i)));
    mat_cell{i} = pauli_mat;
    Var_cell{i} = Var_cell{i}*coeffs_list(i);
end
if nc >0
    H_hk = setup(H_hk,Var_cell,k_cell,mat_cell,silence);
end
end
