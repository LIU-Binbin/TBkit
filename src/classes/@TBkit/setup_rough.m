function H_hk = setup_rough(H_hk,symbolic_polynomial,pauli_mat,silence)
%SETUP_ROUGH Configure HK object from symbolic polynomial and Pauli matrices
%   H_HK = SETUP_ROUGH(H_HK, SYMBOLIC_POLYNOMIAL, PAULI_MAT, SILENCE) sets up
%   Hamiltonian terms by expanding a symbolic polynomial with Pauli matrices.
%
%   Inputs:
%       H_hk              - Input HK object
%       symbolic_polynomial - Symbolic polynomial expression
%       pauli_mat         - Pauli matrix coefficients
%       silence           - Suppress progress display if true (default: false)
%
%   Output:
%       H_hk              - Configured HK object
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
