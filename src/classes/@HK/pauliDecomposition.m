function [H_sym_pauli, H_latex_pauli] = pauliDecomposition(H_hk)
%PAULIDECOMPOSITION Decompose Hamiltonian in Pauli matrix basis
%
% Syntax:
%   [H_sym, H_latex] = pauliDecomposition(H_hk)
%
% Input:
%   H_hk - HK Hamiltonian object
%
% Outputs:
%   H_sym_pauli - Symbolic Pauli matrix decomposition
%   H_latex_pauli - LaTeX representation of decomposition
%
% Description:
%   Performs decomposition of the Hamiltonian into Pauli matrix components:
%   1. Extracts symbolic Hamiltonian (Hk_sym)
%   2. Delegates to TBkit.pauliDecomposition for core operation
%   3. Returns both symbolic and LaTeX formatted results
%
% Note:
%   Requires Hamiltonian to be compatible with Pauli matrix basis
%   Preserves original Hamiltonian structure
%
% Example:
%   [H_pauli, H_tex] = Hk.pauliDecomposition();
H_sym = H_hk.Hk_sym;
[H_sym_pauli,~,H_latex_pauli]= TBkit.pauliDecomposition(H_sym);
end
