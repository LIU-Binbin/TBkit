function [klist_cart,klist_frac,klist_cart_plot,sizemesh,Gk_,Grid] = kmesh2D(Rm,options,optionsEdge)
%KMESH2D Generate 2D k-mesh for surface calculations
%   [k_cart, k_frac, k_plot, sz, Gk, Grid] = kmesh2D(Rm, opts, edgeOpts)
%
%   Inputs:
%       Rm         - Lattice vectors or object containing Rm
%       options    - Mesh parameters:
%           .knum1/2    - Grid dimensions (default: 51)
%           .kstart     - Mesh origin (default: [-0.5,-0.5,0])
%           .kdir1/2    - Mesh directions (default: [1,0,0], [0,1,0])
%           .cartesian  - Use cartesian coordinates (default: false)
%       optionsEdge - Edge handling:
%           .full       - Include edge points (default: false)
%
%   Outputs:
%       klist_cart      - Cartesian k-points
%       klist_frac      - Fractional k-points
%       klist_cart_plot - Shifted points for plotting
%       sizemesh        - Mesh dimensions [N1,N2]
%       Gk_             - Used reciprocal lattice vectors
%       Grid            - 3D mesh grid
arguments
    Rm =POSCAR_read;
    options.knum1   = 51;
    options.knum2   = 51;
    options.kstart  = [-0.5,-0.5,0];
    options.kdir1   = [1,0,0];
    options.kdir2   = [0,1,0];
    options.cartesian = false;
    options.dir_seq = [1,2,3];
    options.dir_start = 'kcar';
    options.shift = true;
    optionsEdge.full = false;
end
if ~isa(Rm,'double') && ~isa(Rm,'sym')
    Rm = Rm.Rm;
end
Gk = (2*pi*eye(size(Rm,1))/Rm).';
if size(Gk,1) < 3 || size(Gk,2) <3
    Gk(3,3) = 1;
end
%
if options.cartesian
    Gk_ = vasplib.CartisianMat(Gk,options.dir_seq,options.dir_start);
    kstart_s  = options.kstart * Gk_ /Gk;
else
    Gk_ = Gk;
    kstart_s = options.kstart;
end
%
knum1 = options.knum1+1;
knum2 = options.knum2+1;
[~,klist_s_1,~,~] =...
    vasplib.kpathgen([[0,0,0];options.kdir1],knum1,Gk_,Gk);
%             klist_l = zeros(size(klist_s_1,1),1);
%             klist_l(1) = sum(sign(klist_s_1(1,:)))*norm(klist_s_1(1,:)*(eye(3)*2*pi));
%klist_s_1_ = klist_s_1 ;
%normklist_l = norm(options.kdir1)/norm(klist_s_1(end,:));
%             for i = 1:size(klist_s_1_,1)
%                 klist_l(i) = norm(klist_s_1_(i,:)*(eye(3)*2*pi))*normklist_l;
%             end
%             klist_l = klist_l + sum(sign(kstart_s))* norm(kstart_s*(eye(3)*2*pi))*normklist_l;
[~,klist_s_2,~,~] =...
    vasplib.kpathgen([[0,0,0];options.kdir2],knum2,Gk_,Gk);
% half-closed half-open
dk_s_1 = options.kdir1/(options.knum1-1);
dk_s_2 = options.kdir2/(options.knum2-1);
dk_s_start = dk_s_1/2 + dk_s_2/2;
if optionsEdge.full
    klist_s_1_ = klist_s_1(1:knum1-1,:);
    klist_s_2_ = klist_s_2(1:knum2-1,:);
    knum1_ = knum1-1;
    knum2_ = knum2-1;
    sizemesh = [knum1_,knum2_];
else
    klist_s_1 = klist_s_1(1:knum1-1,:);
    klist_s_2 = klist_s_2(1:knum2-1,:);
    knum1 = knum1-1;
    knum2 = knum2-1;
    sizemesh = [knum1,knum2];
end
%BC_WAVECAR = zeros(TBkitobj.Basis_num,options.knum1*options.knum2);
%BCCAR = zeros(options.knum1,options.knum_evol);
%BF_WAVECAR = zeros(length(BAND_index),length(BAND_index),options.knum_evol);
%kstart_r = options.kstart*Gk_;
% make klist
%
klist_frac = repmat(klist_s_1,[ knum2 1 ])+kron(klist_s_2,ones( knum1, 1 )) + kstart_s;
klist_cart = klist_frac*Gk_;
if options.shift
    shift = dk_s_start*Gk_;
else
    shift = [0,0,0];
end
if optionsEdge.full
    klist_frac_ = repmat(klist_s_1_,[ knum2_ 1 ])+kron(klist_s_2_,ones( knum1_, 1 )) + kstart_s;
    klist_cart_ = klist_frac_*Gk_;
    klist_cart_plot = klist_cart_ + shift;
else
    klist_cart_plot = klist_cart + shift;
end
Grid(:,:,1)  = reshape(klist_cart_plot(:,1),sizemesh);
Grid(:,:,2)  = reshape(klist_cart_plot(:,2),sizemesh);
Grid(:,:,3)  = reshape(klist_cart_plot(:,3),sizemesh);

end
