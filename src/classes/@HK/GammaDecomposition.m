function [H_sym_Gamma,H_latex_Gamma] = GammaDecomposition(H_hk)
%GAMMADECOMPOSITION Decompose Hamiltonian at Gamma point into symmetry components
%
% Syntax:
%   [H_sym_Gamma,H_latex_Gamma] = GammaDecomposition(H_hk)
%
% Input:
%   H_hk - HK Hamiltonian object
%
% Outputs:
%   H_sym_Gamma - Symbolic decomposition at Γ point (k=0)
%   H_latex_Gamma - LaTeX representation of the decomposition
%
% Description:
%   Performs symmetry decomposition of the Hamiltonian at the Γ point by:
%   1. Extracting the symbolic Hamiltonian (Hk_sym)
%   2. Delegating to TBkit.GammaDecomposition for actual decomposition
%   Returns both symbolic and LaTeX formatted results for analysis/display.
%
% Note:
%   Relies on TBkit's GammaDecomposition for core functionality
%
% Example:
%   [H_sym, H_tex] = Hk_obj.GammaDecomposition();
H_sym = H_hk.Hk_sym;
[H_sym_Gamma,~,H_latex_Gamma] = TBkit.GammaDecomposition(H_sym);
end
