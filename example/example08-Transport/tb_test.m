hr = HR.from_wannier90();
%%
kcube_bulk =  [-0.5 -0.5 -0.5; 1 0 0; 0 1 0; 0 0 1];
[klist_cart, klist_frac] = kcubegen3D('Rm', hr.Rm, 'KCUBE_BULK', kcube_bulk, 'nk', [300 300 1]);

mu_list = linspace(-0.2, 0.2, 201);
chi_mu = SOAHC_int(hr, [1 2 2], klist_frac, mu_list, 'ncore',4 , 'T',20);
%% second-order anomalous Hall effect
if 1 == 0
    chi_mu = chi_mu .* hr.Rm(3,3) .* 0.1; % 3D to 2D, Ang to nm
    
    [fig, ax] = Figs(1,1);
    plot(mu_list, chi_mu, 'LineWidth', 2, 'Color', 'r')
    xlabel("\mu (eV)")
    ylabel("\chi^{int}_{xyy} (nmAV^{-2})")
end