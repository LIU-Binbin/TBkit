function [W,D,dH_dk_xyz] = fft(H_hk, kpoint, precomputed)
arguments
    H_hk HK
    kpoint % cart
    precomputed = [] ; % 可选缓存
end

Hout = H_hk.Hfun(kpoint(1), kpoint(2), kpoint(3));
Hout = (Hout+Hout')/2;
[W,D] = eig(Hout, 'vector');
dH_dk_xyz = H_hk.dH_dk_fun(kpoint(1), kpoint(2), kpoint(3));
end
