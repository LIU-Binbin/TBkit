% HTRIG2HK Convert a Htrig object to a corresponding HK object.
%
% SYNTAX:
%   H_hk = Htrig2HK(H_htrig, kpoints_f, options)
%
% DESCRIPTION:
%   This function transforms a Htrig object (tight-binding representation in real space)
%   to an HK object (tight-binding representation in k-space). The conversion involves a 
%   transformation of the Hamiltonian's symbolic terms (HsymL_trig) into polynomial 
%   approximations (using taylor expansion) centered at shifted k-points computed from the 
%   provided kpoints_f and the reciprocal lattice (Gk). The order of the Taylor expansion 
%   is specified by options.Order.
%
%   If the 'sym' option is true, the reciprocal lattice matrix Rm is converted to symbolic 
%   form and the transformed k-points (kpoints_r) are displayed for verification.
%
% INPUTS:
%   H_htrig    - A Htrig object.
%   kpoints_f  - A 1x3 vector representing the k-point in f-space (default: [0, 0, 0]).
%   options    - A structure with optional fields:
%                   sym   - Logical flag to use symbolic Rm and display kpoints_r (default: false).
%                   Order - Order of the Taylor expansion (default: 1).
%
% OUTPUT:
%   H_hk       - An HK object with Basis_num and Order inherited from H_htrig, containing 
%                the k-space tight-binding representation.
%
% EXAMPLE:
%   % Convert a Htrig object H to its HK representation with default options:
%   H_hk = Htrig2HK(H, [0, 0, 0], struct('sym', false, 'Order', 1));
%
% SEE ALSO:
%   HK, TBkitCopy, taylor, subs, setup_rough, CmdLineProgressBar
%
function H_hk = Htrig2HK(H_htrig, kpoints_f, options)
arguments
    H_htrig Htrig;
    kpoints_f = [0,0,0];
    options.sym = false;
    options.Order = 1;
end

syms k_x k_y k_z real;
H_hk = HK(H_htrig.Basis_num, options.Order);
H_hk = H_hk.TBkitCopy(H_htrig);

if options.sym
    H_hk.Rm = sym(H_hk.Rm);
    kpoints_r = kpoints_f * H_hk.Gk;
    fprintf('please check whether the kpoint (Cartesian) is properly set.\n');
    disp(kpoints_r);
else
    kpoints_r = kpoints_f * H_hk.Gk;
end

pb = TBkit_tool_outer.CmdLineProgressBar('Transforming ');
for i = 1:H_htrig.Kinds
    pb.print(i, H_htrig.Kinds, ' th HsymL_trig into HK obj ...');
    tmp_HsymL_trig = subs(H_htrig.HsymL_trig(i), [k_x k_y k_z], [k_x - kpoints_r(1), k_y - kpoints_r(2), k_z - kpoints_r(3)]);
    symbolic_polynomial = taylor(tmp_HsymL_trig, [k_x k_y k_z], 'Order', options.Order+1);
    H_hk = H_hk.setup_rough(symbolic_polynomial, H_htrig.HcoeL(:,:,i), true);
end
pb.delete();
end
