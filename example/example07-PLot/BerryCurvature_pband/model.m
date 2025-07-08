%% useful_tools
sigma_0 = pauli_matrix(0); sigma_x = pauli_matrix(1); sigma_y = pauli_matrix(2); sigma_z = pauli_matrix(3);
%%
syms k_x k_y k_z real
syms v w d0 d2 real 
%%
H_kp = HK(2,2);
H_kp = H_kp...
     + Term(-v*k_x, sigma_y)...
     + Term( v*k_y, sigma_x)...
     + Term( w*k_x, sigma_0)...
     + Term(-d2*(k_x^2+k_y^2) + d0, sigma_z);
%%
v = 1;
w = 0.5*v;
d0 = 0.1;
d2 = 1;
%%
H_kp_n = H_kp.Subsall();
EIGENCAR = H_kp_n.EIGENCAR_gen();
%%
[klist_cart, klist_l, klist_frac, ~, ~] = kpathgen3D();
nk = size(klist_cart,1);
Omega_xy = zeros(2, nk);
for i = 1:nk
    Omega_xy(:, i) = BerryCurvature_nk(H_kp_n, [1 2], klist_cart(i,:));
end
%%
[~, ax] = Figs(1,1);
WEIGHT = {Omega_xy};
pbandplot(WEIGHT, EIGENCAR, 'Ecut', [-0.5 0.5], 'LineWidth', 4, ...
    'plotMode', 'RWB', 'ax', ax);
colorbar
bandplot(EIGENCAR, [-0.5 0.5], 'LineSpec', '--', 'Color','k', 'ax', ax);