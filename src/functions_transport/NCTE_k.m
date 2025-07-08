function alpha_mu = NCTE_k(Ham, tensor_index, kpoint,mu_list, T, eps)
arguments
    Ham TBkit
    tensor_index (1,3) double
    kpoint double
    mu_list double
    T = 50 % Kelvin
    eps = 1e-4
end
Nbands = Ham.Basis_num;
a = tensor_index(1);
b = tensor_index(2);
c = tensor_index(3);
% T = options.T;
% eps = options.eps;
%%
% [EIG_ki, WAV_ki] = Ham.EIGENCAR_gen('klist', kpoint, 'printmode', false);
% dH_dk_xyz = Ham.dH_dk(kpoint);
[WAV_ki,EIG_ki,dH_dk_xyz] = Ham.fft(kpoint);
 dEnm = EIG_ki - EIG_ki';  % 能级差矩阵
% dEnm = repmat(EIG_ki, 1, Nbands) - repmat(EIG_ki', Nbands, 1);

inv_dEnm_sq = zeros(Nbands);
valid = abs(dEnm) > eps;
% inv_dEnm(~is_degenerated) = 1./dEnm(~is_degenerated);
inv_dEnm_sq(valid) = 1./(dEnm(valid).^2);
%%
VEC_ki = zeros(Nbands, Nbands, 3);
VEC_nn_ki = zeros(Nbands, 3);
for i = 1:3
    VEC_ki(:,:,i) = WAV_ki' * dH_dk_xyz(:,:,i) * WAV_ki;
    VEC_nn_ki(:,i)= real(diag(VEC_ki(:,:,i)));
end
%%
alpha_n = zeros([Nbands,1]);
% % Nbands = 18; vectorized same with JIT
% $ Nbands = 36; JIT better than vectorized
% if b ~= c
%     ImV = imag(VEC_ki(:,:,b) .* conj(VEC_ki(:,:,c)));
%     Omega_n = -2 * sum(ImV .* inv_dEnm_sq, 2);
%     alpha_n = alpha_n + VEC_nn_ki(:,a) .* Omega_n;
% end
% 
% if c ~= a
%     ImV = imag(VEC_ki(:,:,c) .* conj(VEC_ki(:,:,a)));
%     Omega_n = -2 * sum(ImV .* inv_dEnm_sq, 2);
%     alpha_n = alpha_n + VEC_nn_ki(:,b) .* Omega_n;
% end
% 
% if a ~= b
%     ImV = imag(VEC_ki(:,:,a) .* conj(VEC_ki(:,:,b)));
%     Omega_n = -2 * sum(ImV .* inv_dEnm_sq, 2);
%     alpha_n = alpha_n - VEC_nn_ki(:,c) .* Omega_n;
% end



if b ~= c
    Omega_n = zeros([Nbands,1]);
    for n = 1:Nbands
        for m = 1:Nbands
            Omega_n(n) = Omega_n(n) - 2*imag(VEC_ki(n,m,b) * VEC_ki(m,n,c)) * inv_dEnm_sq(n,m);
        end
    end
    alpha_n = alpha_n + VEC_nn_ki(:,a) .* Omega_n;
end

if c ~= a
    Omega_n = zeros([Nbands,1]);
    for n = 1:Nbands
        for m = 1:Nbands
            Omega_n(n) = Omega_n(n) - 2*imag(VEC_ki(n,m,c) * VEC_ki(m,n,a)) * inv_dEnm_sq(n,m);
        end
    end
    alpha_n = alpha_n + VEC_nn_ki(:,b) .* Omega_n;
end

if a ~= b
    Omega_n = zeros([Nbands,1]);
    for n = 1:Nbands
        for m = 1:Nbands
            Omega_n(n) = Omega_n(n) - 2*imag(VEC_ki(n,m,a) * VEC_ki(m,n,b)) * inv_dEnm_sq(n,m);
        end
    end
    alpha_n = alpha_n - VEC_nn_ki(:,c) .* Omega_n;
end
% 费米分布计算
E_minus_mu = EIG_ki - mu_list(:)';  % 自动扩展
f1 = Fermi_1(E_minus_mu, T);
% %%
% Nmu = length(mu_list);
% E_minus_mu = repmat(EIG_ki, 1, Nmu) - repmat(mu_list, Nbands, 1);
% f1 = Fermi_1(E_minus_mu, T);
% 
alpha_mu = tensorprod( alpha_n, f1.*E_minus_mu, 1, 1);
end