function f2 = Fermi_2(E_minus_mu, T)
% second derivative of Fermi-Dirac distribution 
% 二阶导容易出现数值不稳定性，尽量改用差分代替
k_B = 8.6173324e-5; % eV/K
beta = 1/(T * k_B);

f0 = Fermi_0(E_minus_mu, T);
f2  = beta^2 .* f0 .* (1 - f0) .* (1 - 2*f0);
end