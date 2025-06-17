function [kloop1_frac,kloop2_frac,kloop1_cart,kloop2_cart,klist_l,kstart_frac,kstart_cart] = kloop2D(Rm,options)
%KLOOP2D Generate 2D k-loop for integration/evolution calculations
%   [kf1, kf2, kc1, kc2, kl, ksf, ksc] = kloop2D(Rm, opts)
%
%   Inputs:
%       Rm        - Lattice vectors
%       options   - Loop parameters:
%           .knum_evol/int - Number of points
%           .kstart        - Loop center
%           .kevolution    - Evolution direction
%           .kintegral     - Integration direction
%           .cartesian     - Use cartesian coords (default: false)
%
%   Outputs:
%       Fractional and cartesian k-points with length metrics
arguments
    Rm
    options.knum_evol   = 51;
    options.knum_int   = 51;
    options.kstart  = [-0.5,-0.5,0];
    options.kevolution   = [1,0,0];
    options.kintegral   = [0,1,0];
    options.cartesian = false;
    options.dir_seq = [1,2,3];
    options.dir_start = 'kcar';
    options.shift = true;
end
if ~isa(Rm,'double')
    Rm = Rm.Rm;
end
Gk = (2*pi*eye(size(Rm,1))/Rm).';
if size(Gk,1) < 3 || size(Gk,2) <3
    Gk(3,3) = 1;
end
if options.cartesian
    Gk_ = TBkit.CartisianMat(Gk,options.dir_seq,options.dir_start);
    kstart_frac  = options.kstart * Gk_ /Gk;
else
    Gk_ = Gk;
    kstart_frac = options.kstart;
end
%
knum1 = options.knum_evol;
knum2 = options.knum_int;
%
[kloop1_cart,kloop1_frac,~,~] =...
    TBkit.kpathgen([[0,0,0];options.kevolution],knum1,Gk_,Gk);
[kloop2_cart,kloop2_frac,~,~] =...
    TBkit.kpathgen([[0,0,0];options.kintegral],knum2,Gk_,Gk);
%
klist_l = zeros(knum1,1);
normklist_l = norm(options.kevolution)/norm(kloop1_frac(end,:));
for i = 1:knum1
    klist_l(i) = norm(kloop1_frac(i,:)*(eye(3)*2*pi))*normklist_l;
end
klist_l = klist_l + sum(sign(kstart_frac))* norm(kstart_frac*(eye(3)*2*pi))*normklist_l;
kstart_cart = options.kstart*Gk_;
end
