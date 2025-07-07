function H_hk = get_dH_dk_func(H_hk)
syms k_x k_y k_z real;
symL = H_hk.HsymL;

dsym_dkx = diff(symL, k_x);  % Derivative with respect to k_x
dsym_dky = diff(symL, k_y);  % Derivative with respect to k_y
dsym_dkz = diff(symL, k_z);  % Derivative with respect to k_z

nbands = H_hk.Nbands;
HnumL = reshape(H_hk.HnumL, nbands^2, []);

dH_dk_sym = zeros(nbands, nbands, 3);
dH_dk_sym(:,:,1) = reshape(HnumL * dsym_dkx.', nbands, nbands);
dH_dk_sym(:,:,2) = reshape(HnumL * dsym_dky.', nbands, nbands);
dH_dk_sym(:,:,3) = reshape(HnumL * dsym_dkz.', nbands, nbands);

% Convert the symbolic derivatives into MATLAB function handles
H_hk.dH_dk = matlabFunction(dH_dk_sym, 'Vars', [k_x, k_y, k_z]);

% Convert the symbolic derivatives into MATLAB function handles
% HsymL_fun = matlabFunction( H_hk.HsymL,'Vars',H_hk.VarsSeqLcart(1:H_hk.Dim));
H_hk.Hfun = matlabFunction(H_hk.sym, 'Vars', [k_x, k_y, k_z]);
% syms k_x k_y k_z real;
% 
% symL = H_hk.HsymL;
% dsym_dk = zeros(length(symL),3);
% 
% dsym_dk(:,1) = diff(symL, k_x);  % Derivative with respect to k_x
% dsym_dk(:,2) = diff(symL, k_y);  % Derivative with respect to k_y
% dsym_dk(:,3) = diff(symL, k_z);  % Derivative with respect to k_z
% 
% dnum_dk = double(subs(dsym_dk, [k_x, k_y, k_z], kpoint));
% dH_dk_xyz = tensorprod(H_hk.HnumL, dnum_dk, 3, 1);
end