function alpha_mu = NCTE_k(Ham, tensor_index, kpoint, mu_list, options)
arguments
    Ham TBkit
    tensor_index (1,3) double
    kpoint double
    mu_list double
    options.T = 50 % Kelvin
    options.eps = 1e-4
end
Nbands = Ham.Nbands;
a = tensor_index(1);
b = tensor_index(2);
c = tensor_index(3);
%%
[EIG_ki, WAV_ki] = Ham.EIGENCAR_gen('klist', kpoint);

dEnm = repmat(EIG_ki, 1, Nbands) - repmat(EIG_ki', Nbands, 1);
inv_dEnm = zeros(Nbands, Nbands);
is_degenerated = abs(dEnm) < options.eps;
inv_dEnm(~is_degenerated) = 1./dEnm(~is_degenerated);
%%
dH_dk_xyz = Ham.dH_dk(kpoint);
VEC_ki = zeros(Nbands, Nbands, 3);
VEC_nn_ki = zeros(Nbands, 3);
for i = 1:3
    VEC_ki(:,:,i) = WAV_ki' * dH_dk_xyz(:,:,i) * WAV_ki;
    VEC_nn_ki(:,i)= real(diag(VEC_ki(:,:,i)));
end
%%
alpha_n = zeros([Nbands,1]);

if b ~= c
    Omega_n = zeros([Nbands,1]);
    for n = 1:Nbands
        for m = 1:Nbands
            Omega_n(n) = Omega_n(n) - 2*imag(VEC_ki(n,m,b) * VEC_ki(m,n,c)) * inv_dEnm(n,m)^2;
        end
    end
    alpha_n = alpha_n + VEC_nn_ki(:,a) .* Omega_n;
end

if c ~= a
    Omega_n = zeros([Nbands,1]);
    for n = 1:Nbands
        for m = 1:Nbands
            Omega_n(n) = Omega_n(n) - 2*imag(VEC_ki(n,m,c) * VEC_ki(m,n,a)) * inv_dEnm(n,m)^2;
        end
    end
    alpha_n = alpha_n + VEC_nn_ki(:,b) .* Omega_n;
end

if a ~= b
    Omega_n = zeros([Nbands,1]);
    for n = 1:Nbands
        for m = 1:Nbands
            Omega_n(n) = Omega_n(n) - 2*imag(VEC_ki(n,m,a) * VEC_ki(m,n,b)) * inv_dEnm(n,m)^2;
        end
    end
    alpha_n = alpha_n - VEC_nn_ki(:,c) .* Omega_n;
end
%%
Nmu = length(mu_list);
E_minus_mu = repmat(EIG_ki, 1, Nmu) - repmat(mu_list, Nbands, 1);
f1 = Fermi_1(E_minus_mu, options.T);

alpha_mu = tensorprod( alpha_n, f1.*E_minus_mu, 1, 1);
end