function [BCCAR,Grid,klist_r_plot] = BerryCuvature_fun(Hfun,optionsK,options,optionsPlot)
%BERRYCUVATURE_FUN Calculate Berry curvature from function handle
%
%   Syntax:
%       [BCCAR,Grid,klist_r_plot] = BerryCuvature_fun(Hfun,optionsK,options,optionsPlot)
%
%   Description:
%       Computes Berry curvature distribution over a 2D k-mesh from a
%       Hamiltonian function handle. Supports both numerical and symbolic
%       inputs.
%
%   Inputs:
%       Hfun      - Hamiltonian function handle or symbolic expression
%       optionsK  - k-mesh options (same as BC_2D)
%       options   - Calculation options:
%                   knodes     - Number of k-points
%                   Rm         - Real-space lattice vectors
%                   plot       - Plot flag
%                   car        - Cartesian coordinates flag
%                   BAND_index - Band indices
%                   BCfun      - Direct Berry curvature function flag
%       optionsPlot - Plotting options
%
%   Outputs:
%       BCCAR      - Berry curvature array
%       Grid       - k-mesh grid
%       klist_r_plot - k-points for plotting
arguments
    Hfun ;
    optionsK.knum1   = 51;
    optionsK.knum2   = 51;
    optionsK.kstart  = [-0.5,-0.5,0];
    optionsK.kdir1   = [1,0,0];
    optionsK.kdir2   = [0,1,0];
    optionsK.cartesian = false;
    optionsK.dir_seq = [1,2,3];
    optionsK.dir_start = 'kcar';
    options.knodes = 100;
    options.Rm = [];
    options.plot = false;
    options.car = true;
    options.BAND_index = [];
    options.BCfun = false;
    optionsPlot.oneshot = true;
    optionsPlot.view = [0,90];
end
if isempty(options.Rm)
    Rm = POSCAR_read;
else
    Rm = options.Rm;
end
Gk = (2*pi*eye(3)/Rm).';
%
optionsKcell = namedargs2cell(optionsK);
[klist_cart,~,klist_r_plot,sizemesh,Gk_,Grid] = TBkit.kmesh2D(Rm,optionsKcell{:},'full',true);
if isa(Hfun,'sym')
    Hfun = matlabFunction(Hfun,'Vars',[sym('k_x'),sym('k_y'),sym('k_z')]);
end
if ~options.BCfun
    testH = Hfun(0,0,0);
    Norb = size(testH,1);
    if isempty(options.BAND_index)
        set_divide = 2;
        BAND_index = 1:(Norb/set_divide);
    else
        BAND_index = options.BAND_index;
    end
    % integral method
    [~,BC_WAVECAR] = TBkit.EIGENSOLVE(Hfun,klist_cart,Norb);
    BC_WAVECAR = BC_WAVECAR(:,BAND_index,:);
    %reshape(BC_WAVECAR(:,BAND_index,:),[TBkitobj.Basis_num options.knum1*options.knum2]);
    BCCAR = TBkit.BerryCuvature_2D( BC_WAVECAR ,sizemesh);
else
    BCfun = Hfun;
    BCCAR = zeros(sizemesh);
    for ki = 1:numel(BCCAR)
        K = klist_r_plot(ki,:);
        Bcki = BCfun(K(1),K(2),K(3));
        BCCAR(ki) = sum(Bcki);
    end
end
if options.plot
    knum1 = sizemesh(1);
    knum2 = sizemesh(2);
    dk_1 = (optionsK.kdir1)/knum1*Gk_;
    dk_2 = (optionsK.kdir2)/knum2*Gk_;
    if ~options.BCfun
        dS = 1;
    else
        dS = cross(dk_1,dk_2);
    end
    optionsPlotcell = namedargs2cell(optionsPlot);
    [~,~] = TBkit.ShowSurfIntegral(Gk,klist_cart,dk_1,dk_2,BCCAR(:)*norm(dS)/(2*pi),optionsPlotcell{:});
end
end