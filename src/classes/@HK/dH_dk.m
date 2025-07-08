function dH_dk_xyz = dH_dk(H_hk, kpoint)
dH_dk_xyz = H_hk.dH_dk_fun(kpoint(1), kpoint(2), kpoint(3));
end