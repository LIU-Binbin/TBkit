function Omega_ab = BerryCurvature_nk(Ham, tensor_index, kpoint, options)
arguments
    Ham TBkit
    tensor_index (1,2) double
    kpoint (1,3) double
    options.eps = 1e-4
end
Nbands = Ham.Basis_num;
a = tensor_index(1);
b = tensor_index(2);
%%
[EIG_ki, WAV_ki] = Ham.EIGENCAR_gen('klist', kpoint);

dEnm = repmat(EIG_ki, 1, Nbands) - repmat(EIG_ki', Nbands, 1);
inv_dEnm = zeros(Nbands, Nbands);
is_degenerated = abs(dEnm) < options.eps;
inv_dEnm(~is_degenerated) = 1./dEnm(~is_degenerated);
%%
dH_dk_xyz = Ham.dH_dk(kpoint);
VEC_ki = zeros(Nbands, Nbands, 3);
for i = 1:3
    VEC_ki(:,:,i) = WAV_ki' * dH_dk_xyz(:,:,i) * WAV_ki;
end
%%
Omega_ab = zeros([Nbands,1]);

for n = 1:Nbands
    for m = 1:Nbands
        Omega_ab(n) = Omega_ab(n) - 2*imag(VEC_ki(n,m,a) * VEC_ki(m,n,b)) * inv_dEnm(n,m)^2;
    end
end
end