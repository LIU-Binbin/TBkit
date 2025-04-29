function [BFCAR,BF_WAVECAR,klist_l] = WilsonLoop_fun(Hfun,optionsK,options)
% WILSONLOOP_FUN Calculate Wilson loop spectrum for a given Hamiltonian
%
% Syntax:
%   [BFCAR,BF_WAVECAR,klist] = WilsonLoop_fun(Hfun,optionsK)
%   [BFCAR,BF_WAVECAR,klist] = WilsonLoop_fun(...,options)
%
% Description:
%   Computes Wilson loop spectrum (Berry phases) for a 2D k-loop
%   parameterization of a given Hamiltonian function.
%
% Input Arguments:
%   Hfun - Hamiltonian function handle or symbolic matrix
%   optionsK - k-path parameters (same as WannierCenter)
%   options - Additional options (same as WannierCenter)
%
% Output Arguments:
%   BFCAR - Berry phase spectrum
%   BF_WAVECAR - Corresponding eigenvectors
%   klist_l - k-point list used
%
% Example:
%   H = @(kx,ky) [0,exp(1i*kx)+exp(1i*ky); exp(-1i*kx)+exp(-1i*ky),0];
%   [bf,vec] = WilsonLoop_fun(H,struct('knum_int',50));
arguments
    Hfun ;
    optionsK.knum_int    = 31;
    optionsK.knum_evol   = 51;
    optionsK.kstart      = [0,0,0];
    optionsK.kintegral   = [0,1,0];
    optionsK.kevolution  = [1,0,0];
    optionsK.cartesian = false;
    optionsK.dir_seq = [1,2,3];
    optionsK.dir_start = 'kcar';
    options.WLfun = false;
    options.ax = handle([]);
    options.plot = false;
    options.V = [];
    options.Accuracy = 1E-12;
    options.LWAVE = false;
    options.BAND_index = [];
    options.script {mustBeMember(options.script,...
        {...
        'nu_x(k_y)','nu_y(k_x)',....
        'nu_y(k_z)','nu_z(k_y)',...
        'nu_z(k_x)','nu_x(k_z)',...
        'nu_1(k_2)','nu_2(k_1)',...
        'nu_2(k_3)','nu_3(k_2)',...
        'nu_3(k_1)','nu_1(k_3)',...
        ''})}= '';
end
% --- nargin
switch class(Hfun)
    case {'sym'}
        syms k_x k_y k_z real;
        Hfun = matlabFunction(Hfun,'Vars',[k_x,k_y,k_z]);
    case 'handle'

    otherwise
end
if ~options.WLfun
    testH = Hfun(0,0,0);
    NBAND = size(testH,1);
    if isempty(options.BAND_index)
        set_divide = 2;
        BAND_index = 1:(Norb/set_divide);
    else
        BAND_index = options.BAND_index;
    end
else
end
% use script
switch options.script
    case {'nu_x(k_y)','nu_1(k_2)'}
        optionsK.kintegral   = [1,0,0];
        optionsK.kevolution  = [0,1,0];
        optionsK.kstart      = [0,-0.5,0];
    case {'nu_y(k_z)','nu_2(k_3)'}
        optionsK.kintegral   = [0,1,0];
        optionsK.kevolution  = [0,0,1];
        optionsK.kstart      = [0,0,-0.5];
    case {'nu_z(k_x)','nu_3(k_1)'}
        optionsK.kintegral   = [0,0,1];
        optionsK.kevolution  = [1,0,0];
        optionsK.kstart      = [-0.5,0,0];
    case {'nu_y(k_x)','nu_2(k_1)'}
        optionsK.kintegral   = [0,1,0];
        optionsK.kevolution  = [1,0,0];
        optionsK.kstart      = [-0.5,0,0];
    case {'nu_z(k_y)','nu_3(k_2)'}
        optionsK.kintegral   = [0,0,1];
        optionsK.kevolution  = [0,1,0];
        optionsK.kstart      = [0,-0.5,0];
    case {'nu_x(k_z)','nu_1(k_3)'}
        optionsK.kintegral   = [1,0,0];
        optionsK.kevolution  = [0,0,1];
        optionsK.kstart      = [0,0,-0.5];
    case {''}
end
%
if contains(string(options.script),["x","y","z"])
    optionsK.cartesian = true;
    optionsK.dir_seq = [1,2,3];
    optionsK.dir_start = 'kcar';
end
optionsKcell = namedargs2cell(optionsK);
% Gen kloop 2D
[~,~,kloop1_cart,kloop2_cart,klist_l,~,kstart_cart] = TBkit.kloop2D(TBkitobj.Rm,optionsKcell{:});
% prepare
BFCAR = zeros(length(BAND_index),optionsK.knum_evol);
BF_WAVECAR = zeros(length(BAND_index),length(BAND_index),optionsK.knum_evol);
if options.LWAVE
    if isempty(options.V )
        V = diag((1:NBAND) *sqrt(options.Accuracy));
    else
        V = options.V ;
    end
    %EIGENCAR = zeros(length(BAND_index),optionsK.knum_int,optionsK.knum_evol);
    WAVELOOPCAR = zeros(NBAND,length(BAND_index),optionsK.knum_int,optionsK.knum_evol);
end
for i = 1:optionsK.knum_evol
    klist_tmp = kloop1_cart(i,:)+kloop2_cart+kstart_cart;
    [E,WAVECAR_loop] = TBkit.EIGENSOLVE(Hfun, ...
        klist_tmp,NBAND);
    WAVECAR_loop_tmp = WAVECAR_loop(:,BAND_index,:);
    % The last bloch state is the same as the first up to a phase factor
    WAVECAR_loop_tmp(:,:,end) = WAVECAR_loop_tmp(:,:,1);
    if options.LWAVE
        WAVECAR_loop_tmp = TBkit.cleanWAVECAR(WAVECAR_loop_tmp,E(BAND_index,:),V,options.Accuracy);
        WAVELOOPCAR(:,:,:,i) = WAVECAR_loop_tmp;
    end
    [BFCAR(:,i),BF_WAVECAR(:,:,i)] = TBkit.wancenter_1D(WAVECAR_loop_tmp);
end

end