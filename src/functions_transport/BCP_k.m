function Gbc = BCP_k(Ham, tensor_index, kpoint,selectbands)
arguments
    Ham TBkit
    tensor_index (1,2) double
    kpoint double
    selectbands
end
Nbands = Ham.Basis_num;
b = tensor_index(1);
c = tensor_index(2);
% T = options.T;
% eps = options.eps;
%%
[~,EIG_ki,dH_dk_xyz] = Ham.fft(kpoint);
% dH_dk_xyz = Ham.dH_dk(kpoint);
 dEnm = EIG_ki - EIG_ki';  % 能级差矩阵
% dEnm = repmat(EIG_ki, 1, Nbands) - repmat(EIG_ki', Nbands, 1);

inv_dEnm_tr = zeros(Nbands);
valid = abs(dEnm) > eps;
% inv_dEnm(~is_degenerated) = 1./dEnm(~is_degenerated);
inv_dEnm_tr(valid) = 2./(dEnm(valid).^3);
%%
Vb = dH_dk_xyz(:,:,b);
Vc = dH_dk_xyz(:,:,c);
% $G_{b c, n}=2 \operatorname{Re} \sum_{\ell \neq n} \frac{\left\langle v_b\right\rangle_{n \ell}\left\langle v_c\right\rangle_{\ell n}}{\left(\varepsilon_n-\varepsilon_{\ell}\right)^3}$
Gbc = zeros(1,length(selectbands));
for i = 1:length(selectbands)
    n = selectbands(i);
    for l = 1:Nbands
        if l ~= n && inv_dEnm_tr(l,n) ~=0
            Gbc(i) = Gbc(i) + real(Vb(n,l)*Vc(l,n)* inv_dEnm_tr(n,l));
        end
        
    end
end

end