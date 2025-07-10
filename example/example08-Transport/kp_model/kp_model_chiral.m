%% useful_tools
sigma_0 = pauli_matrix(0); sigma_x = pauli_matrix(1); sigma_y = pauli_matrix(2); sigma_z = pauli_matrix(3);
tau_0   = pauli_matrix(0); tau_x   = pauli_matrix(1); tau_y   = pauli_matrix(2); tau_z   = pauli_matrix(3);
%%
syms k_x k_y k_z real
syms v1 v2 v3 w real 
%%
H_kp = HK(2,2);
H_kp = H_kp...
     + Term(v1*k_x, sigma_x)...
     + Term(v2*k_y, sigma_y)...
     + Term(v3*k_z, sigma_z)...
     + Term(w*k_x^2, sigma_0)...
     + Term(w*k_y^2, sigma_0)...
     + Term(w*k_z^2, sigma_0);
%%
v1 = 1;
v2 = 1;
v3 = 1;
delta = 40e-3;
alpha = 1.2;%m* = 1.2 m_e, alpha = 1.2
w = 3.81/alpha;% % $E(k)[\mathrm{eV}] \approx \frac{3.81}{\alpha} k^2$​
% kpoints 
% Gk = 2pi*10();1/nm-1
% frac = -4/(2pi*10)
%%
H_kp_n = H_kp.Subsall();
H_kp_n.bandplot([-0.2 0.6]);

%%
krange = 0.06; % just around the Gamma point
kcube_bulk = krange .* [-0.5 -0.5 -0.5; 1 0 0; 0 1 0; 0 0 1];
[klist_cart, klist_frac] = kcubegen3D('Rm', H_kp_n.Rm, 'KCUBE_BULK', kcube_bulk, 'nk', [300 300 300]);

mu_list = linspace(-0.02, 0.02, 201);
%%
alpha_mu = NCTE(H_kp_n, [3 1 2], klist_cart, mu_list, 'ncore',8 , 'T',20,'batch_size',1e6);
% ~ 6min
kcube_ratio = krange^3;
alpha_mu = alpha_mu .* kcube_ratio;

[fig, ax] = Figs(1,1);
plot(mu_list, alpha_mu, 'LineWidth', 2, 'Color', 'r')
xlabel("\mu (eV)")
ylabel("\alpha/\tau (AV^{-1}K^{-1}s^{-1})")
