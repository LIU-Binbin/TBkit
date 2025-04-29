function [BFCAR,BF_WAVECAR,klist_l,WAVELOOPCAR] = WilsonLoop(TBkitobj,optionsK,options)
% WILSONLOOP Calculate Wilson loop spectrum for tight-binding objects
%
% Syntax:
%   [BFCAR,BF_WAVECAR,klist,WAVELOOP] = WilsonLoop(TBkitobj,optionsK)
%   [BFCAR,BF_WAVECAR,klist,WAVELOOP] = WilsonLoop(...,options)
%
% Description:
%   Computes Wilson loop spectrum for tight-binding objects (Htrig, HK, HR),
%   providing topological information through Berry phase calculations along
%   a 2D k-loop parameterization.
%
% Input Arguments:
%   TBkitobj - Tight-binding object (Htrig, HK, or HR type)
%   optionsK - Structure containing k-path parameters (same as WilsonLoop_fun)
%   options - Additional options (same as WilsonLoop_fun)
%
% Output Arguments:
%   BFCAR - Berry phase spectrum (Nbands x Nk_evol)
%   BF_WAVECAR - Corresponding eigenvectors (Nbands x Nbands x Nk_evol)
%   klist_l - k-point list used in calculation
%   WAVELOOPCAR - Wavefunction data if LWAVE=true (optional)
%
% Example:
%   [wcc,vec,k] = WilsonLoop(Hk,struct('knum_int',50));
%
% See also: WilsonLoop_fun, WannierCenter
arguments
    TBkitobj ;
    optionsK.knum_int    = 31;
    optionsK.knum_evol   = 51;
    optionsK.kstart      = [0,0,0];
    optionsK.kintegral   = [0,1,0];
    optionsK.kevolution  = [1,0,0];
    optionsK.cartesian = false;
    optionsK.dir_seq = [1,2,3];
    optionsK.dir_start = 'kcar';
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
switch class(TBkitobj)
    case {'Htrig','HK'}
        NBAND = TBkitobj.Basis_num;
    case 'HR'
        NBAND = TBkitobj.WAN_NUM;
end
if isempty(options.BAND_index)
    BAND_index = 1:(NBAND/2);
else
    BAND_index = options.BAND_index;
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
[kloop1_frac,kloop2_frac,kloop1_cart,kloop2_cart,klist_l,kstart_frac,kstart_cart] = TBkit.kloop2D(TBkitobj.Rm,optionsKcell{:});
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
switch class(TBkitobj)
    case {'Htrig','HK'}
        if isempty(TBkitobj.orbL)
            TBkitobj.orbL = zeros(TBkitobj.Basis_num,3);
        end
        for i = 1:optionsK.knum_evol
            klist_tmp = kloop1_cart(i,:)+kloop2_cart+kstart_cart;
            [E,WAVECAR_loop] = TBkitobj.EIGENCAR_gen(...
                'klist',klist_tmp,'printmode',false);
            WAVECAR_loop_tmp = WAVECAR_loop(:,BAND_index,:);
            % The last bloch state is the same as the first up to a phase factor
            WAVECAR_loop_tmp(:,:,end) = WAVECAR_loop_tmp(:,:,1).* exp(-1i*(TBkitobj.orbL*TBkitobj.Rm*(klist_tmp(end,:)-klist_tmp(1,:)).'));
            if options.LWAVE
                WAVECAR_loop_tmp = TBkit.cleanWAVECAR(WAVECAR_loop_tmp,E(BAND_index,:),V,options.Accuracy);
                WAVELOOPCAR(:,:,:,i) = WAVECAR_loop_tmp;
            end
            [BFCAR(:,i),BF_WAVECAR(:,:,i)] = TBkit.wancenter_1D(WAVECAR_loop_tmp);
        end
    case 'HR'
        for i = 1:optionsK.knum_evol
            klist_tmp = kloop1_frac(i,:)+kloop2_frac+kstart_frac;
            [E,WAVECAR_loop] = TBkitobj.EIGENCAR_gen(...
                'klist',klist_tmp,...
                'convention','II','printmode',false);
            % If we use convention II, each wavefactor should
            % give back the factor
            % C^{nk}_j = C^{nk}_j_tilde * e^{-ik·tj}.
            WAVECAR_loop_tmp = WAVECAR_loop(:,BAND_index,:);
            % normalize phases to get u instead of phi
            for j =1:size(WAVECAR_loop_tmp,3)
                WAVECAR_loop_tmp(:,:,j) = WAVECAR_loop_tmp(:,:,j).* exp(-2i*pi*(TBkitobj.orbL*klist_tmp(j,:).'));
            end
            % The last bloch state is the same as the first up to a phase factor
            WAVECAR_loop_tmp(:,:,end) = WAVECAR_loop_tmp(:,:,1).* exp(-2i*pi*(TBkitobj.orbL*(klist_tmp(end,:)-klist_tmp(1,:)).'));
            if options.LWAVE
                WAVECAR_loop_tmp = TBkit.cleanWAVECAR(WAVECAR_loop_tmp,E(BAND_index,:),V,options.Accuracy);
                WAVELOOPCAR(:,:,:,i) = WAVECAR_loop_tmp;
            end
            [BFCAR(:,i),BF_WAVECAR(:,:,i)] = TBkit.wancenter_1D(WAVECAR_loop_tmp);
        end
end
end