function [H_sym_pauli,H_latex_pauli] = pauliDecomposition(H_hk)
H_sym = H_hk.Hk_sym;
[H_sym_pauli,~,H_latex_pauli]= TBkit.pauliDecomposition(H_sym);
end
