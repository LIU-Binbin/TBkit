function [W,D,dH_dk_xyz,dH_dk_dk_xyz] = fft_2(H_hk, kpoint, precomputed)
arguments
    H_hk HK
    kpoint % cart
    precomputed = [] ; % 可选缓存
end
Hout = H_hk.Hfun(kpoint(1), kpoint(2), kpoint(3));
Hout = (Hout+Hout')/2;
Hk = (Hout + Hout') / 2;  % Hermitian 修正
[W,D] = eig(Hk, 'vector');
dH_dk_xyz = H_hk.dH_dk_fun(kpoint(1), kpoint(2), kpoint(3));
    % for j = 1:3
    %     dH_dk_dk_R(:,:,i,j) = 1i * tensorprod(HnumList , phase.*vectorList_R(:,i).*vectorList_R(:,j), 3, 1);
    % end
end