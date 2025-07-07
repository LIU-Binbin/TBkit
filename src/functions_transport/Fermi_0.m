function f0 = Fermi_0(E_minus_mu, T)
% Fermi-Dirac distribution
k_B = 8.6173324e-5; % eV/K
beta = 1/(T * k_B);

exp_power = E_minus_mu * beta;

mask = exp_power < 50;

f0 = zeros(size(E_minus_mu));
f0(mask) = 1./(exp(exp_power(mask)) + 1);
end