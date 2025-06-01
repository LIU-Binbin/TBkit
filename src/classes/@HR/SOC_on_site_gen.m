function [H_soc_sym, lambda_syms] = SOC_on_site_gen(H_hr)
[H_soc_sym, lambda_syms] = SOC_on_site(H_hr.orbL, H_hr.elementL, H_hr.quantumL, "mode","basis");
end