Ef = 4.0819;
EIGENCAR = EIGENVAL_read() - Ef;

[WEIGHTCAR_Gd, ~, kpath] = WEIGHTCAR_read("PBAND_Gd_SOC.dat");
WEIGHTCAR_Gd_sum = sum(WEIGHTCAR_Gd,3)/size(WEIGHTCAR_Gd,3);
clear WEIGHTCAR_Gd

[WEIGHTCAR_P, ~, ~] = WEIGHTCAR_read("PBAND_P_SOC.dat");
WEIGHTCAR_P_sum = sum(WEIGHTCAR_P,3)/size(WEIGHTCAR_P,3);
clear WEIGHTCAR_P

[WEIGHTCAR_Zn, ~, ~] = WEIGHTCAR_read("PBAND_Zn_SOC.dat");
WEIGHTCAR_Zn_sum = sum(WEIGHTCAR_Zn,3)/size(WEIGHTCAR_Zn,3);
clear WEIGHTCAR_Zn
%%
WEIGHTCAR = {WEIGHTCAR_Zn_sum, WEIGHTCAR_P_sum, WEIGHTCAR_Gd_sum}; % R G B

pbandplot(WEIGHTCAR, EIGENCAR, 'Ecut', [-2 2], 'LineWidth', 1.5, ...
    'plotMode', 'RGB', ...
    'legend', ["Zn", "P", "Gd"]);

