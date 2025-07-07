function f1 = Fermi_1(E_minus_mu, T)
% first derivative of Fermi-Dirac distribution 
k_B = 8.6173324e-5; % eV/K
beta = 1/(T * k_B);

f0 = Fermi_0(E_minus_mu, T);
f1  = -beta .* f0 .* (1 - f0);
end