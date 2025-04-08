function [H_sym_Gamma,H_latex_Gamma] = GammaDecomposition(H_hk)
H_sym = H_hk.Hk_sym;
[H_sym_Gamma,~,H_latex_Gamma] = TBkit.GammaDecomposition(H_sym);
end
