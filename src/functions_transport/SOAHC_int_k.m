function chi_abc_mu = SOAHC_int_k(Ham, tensor_index, kpoint, mu_list, options)
arguments
    Ham TBkit
    tensor_index (1,3) double
    kpoint (1,3) double
    mu_list double
    options.T = 50 % Kelvin
    options.eps = 1e-4
end
Nbands = Ham.Basis_num; % do not call Ham.Nbands, it has if so it is slower
a = tensor_index(1);
b = tensor_index(2);
c = tensor_index(3);
%%
% H = Ham.Hfun(kpoint(1),kpoint(2),kpoint(3));
% [WAV_ki,EIG_ki_ ]  = eig((H+H')/2);
% EIG_ki = diag(EIG_ki_);
[EIG_ki, WAV_ki] = Ham.EIGENCAR_gen('klist', kpoint, 'printmode', false);

dEnm = repmat(EIG_ki, 1, Nbands) - repmat(EIG_ki.', Nbands, 1);
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
G_ac = zeros([Nbands,1]);
G_bc = zeros([Nbands,1]);

for n = 1:Nbands
    for m = 1:Nbands
        G_ac(n) = G_ac(n) + 2*real(VEC_ki(n,m,a) * VEC_ki(m,n,c)) * inv_dEnm(n,m)^3;
        G_bc(n) = G_bc(n) + 2*real(VEC_ki(n,m,b) * VEC_ki(m,n,c)) * inv_dEnm(n,m)^3;
    end
end

chi_abc = G_bc.*VEC_nn_ki(:,a) - G_ac.*VEC_nn_ki(:,b);

Nmu = length(mu_list);
E_minus_mu = repmat(EIG_ki, 1, Nmu) - repmat(mu_list, Nbands, 1);
f1 = Fermi_1(E_minus_mu, options.T);

chi_abc_mu = tensorprod(chi_abc, f1, 1, 1);
end