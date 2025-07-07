%% useful_tools
sigma_0 = pauli_matrix(0); sigma_x = pauli_matrix(1); sigma_y = pauli_matrix(2); sigma_z = pauli_matrix(3);
tau_0   = pauli_matrix(0); tau_x   = pauli_matrix(1); tau_y   = pauli_matrix(2); tau_z   = pauli_matrix(3);
%%
syms k_x k_y k_z real
syms v1 v2 v3 w delta real 
%% PT-symmetric model
H_kp = HK(4,1);
H_kp = H_kp...
     + Term( w*k_x, tau_0*sigma_0)...
     + Term(v1*k_x, tau_x*sigma_0)...
     + Term(v2*k_y, tau_y*sigma_x)...
     + Term(delta , tau_z*sigma_0);
%%
v1 = 1e6 *constants.hbar_eV_s*1e10;
v2 = 1e6 *constants.hbar_eV_s*1e10;
v3 = 0;
delta = 40e-3;
w = 0.4*v1;
%%
H_kp_n = H_kp.Subsall();
% H_kp_n.bandplot([-0.2 0.2]);
%%
krange = 0.02; % just around the Gamma point
kcube_bulk = krange .* [-0.5 -0.5 -0.5; 1 0 0; 0 1 0; 0 0 1];
[klist_cart, klist_frac] = kcubegen3D('Rm', H_kp_n.Rm, 'KCUBE_BULK', kcube_bulk, 'nk', [300 300 1]);

mu_list = linspace(-0.2, 0.2, 201);
%% second-order anomalous Hall effect
if 1 == 1
    chi_mu = SOAHC_int(H_kp_n, [1 2 2], klist_cart, mu_list, 'ncore',4 , 'T',20);
    
    kcube_ratio = krange^2;
    chi_mu = chi_mu .* kcube_ratio;
    chi_mu = chi_mu .* H_kp_n.Rm(3,3) .* 0.1; % 3D to 2D, Ang to nm
    
    [fig, ax] = Figs(1,1);
    plot(mu_list, chi_mu, 'LineWidth', 2, 'Color', 'r')
    xlabel("\mu (eV)")
    ylabel("\chi^{int}_{xyy} (nmAV^{-2})")
end
%%
