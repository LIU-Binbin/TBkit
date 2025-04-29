function [nested_BFCAR,nested_BF_ALL,klist_l] = nested_WilsonLoop(TBkitobj,optionsK,options,optionsNested)
%NESTED_WILSONLOOP Compute nested Wilson loop spectrum
%
%   Syntax:
%       [nested_BFCAR,nested_BF_ALL,klist_l] = nested_WilsonLoop(TBkitobj,optionsK,options,optionsNested)
%
%   Description:
%       Calculates nested Wilson loop eigenvalues for higher-order topology
%       analysis. Implements the "Wilson loop of Wilson loops" approach.
%
%   Inputs:
%       TBkitobj    - Tight-binding model object
%       optionsK    - k-mesh options for inner loop
%       options     - Calculation options
%       optionsNested - Options for nested loop
%
%   Outputs:
%       nested_BFCAR  - Nested Wilson loop eigenvalues
%       nested_BF_ALL - Full nested Wilson loop data
%       klist_l       - k-point path coordinates
%
%   References:
%       Benalcazar et al., Science 357.6346 (2017): 61-66
%
%   See also: WilsonLoop, WannierCenter
arguments
    TBkitobj ;
    optionsK.knum_int    = 31;
    optionsK.kintegral   = [0,1,0];
    optionsK.kevolution  = [1,0,0];
    optionsK.cartesian = false;
    optionsK.dir_seq = [1,2,3];
    optionsK.dir_start = 'kcar';
    options.BAND_index = [];
    options.ax = handle([]);
    options.plot = false;
    options.V = [];
    optionsNested.kstart      = [0,0,0];
    optionsNested.knum_evol   = 51;
    optionsNested.nested_kevolution  = [0,0,1];
    optionsNested.nested_BAND_index = [];
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
if isempty(optionsNested.nested_BAND_index)
    nested_BAND_index = 1:(NBAND/4);
else
    nested_BAND_index = options.nested_BAND_index;
end
%
%kstart_r = options.kstart*Gk_;
%
optionscell = namedargs2cell(options);
optionsK.knum_evol = optionsK.knum_int;
kstart_s = optionsNested.kstart;
knum_evol = optionsNested.knum_evol;
% nested-1 fix kz semi-fix kx integral ky
% k_z
[~,klist_s_3,~,~] =...
    TBkit.kpathgen([[0,0,0];optionsNested.nested_kevolution],knum_evol,TBkitobj.Gk,TBkitobj.Gk);
klist_l = zeros(size(klist_s_3,1),1);
normklist_l = norm(optionsNested.nested_kevolution)/norm(klist_s_3(end,:));
for i = 1:size(klist_s_3,1)
    klist_l(i) = norm(klist_s_3(i,:)*(eye(3)*2*pi))*normklist_l;
end
klist_l = klist_l + sum(sign(kstart_s))* norm(kstart_s*(eye(3)*2*pi))*normklist_l;
klist_s_3 = klist_s_3+kstart_s;
%
nested_BFCAR = zeros(length(nested_BAND_index),knum_evol);
nested_BF_ALL = zeros(length(nested_BAND_index),optionsK.knum_int,knum_evol);
%
pb = TBkit_tool_outer.CmdLineProgressBar('Nested Berry Phase caculating ');
for k = 1:knum_evol
    % first WAN(y)
    % 在固定kx的情况下，计算沿着ky方向的Wilson loop并得到对应的本征矢量
    % 将先沿着ky方向再沿着kx方向的所有哈密顿量本征波函数存储
    % 注意ky首尾相连 而kx没有
    optionsKbk = optionsK;
    optionsKbk.kstart = klist_s_3(k,:);
    optionsKbk.knum_evol = optionsK.knum_int;
    optionsKbkcell = namedargs2cell(optionsKbk);
    [~,BF_WAVECAR,~,WAVELOOPCAR] = TBkitobj.WilsonLoop(optionsKbkcell{:},optionscell{:},'V',options.V,'LWAVE',true,'LWAVE',true);
    %OccupyBand = size(BFCAR,1);
    NBAND = size(WAVELOOPCAR,1);
    knum_evol_nested = optionsK.knum_int;
    knum_int_nested  = optionsK.knum_int;
    pmulist = zeros(length(nested_BAND_index),knum_evol_nested);
    for ik_evolution = 1:knum_evol_nested %(k_y direction)
        NuWAVELOOP = zeros(NBAND,length(nested_BAND_index),optionsK.knum_int);%(k_x direction)
        for ik_int= 1:knum_int_nested %(k_x direction)
            WAVELOOP_ike_ikint = WAVELOOPCAR(:,:,ik_evolution,ik_int);
            WCCWAVECAR_ikint = BF_WAVECAR(:,nested_BAND_index,ik_int);
            NuWAVELOOP(:,:,ik_int)=sum(WAVELOOP_ike_ikint * WCCWAVECAR_ikint,2);
        end
        % The last  is the same as the first
        NuWAVELOOP(:,:,end) = NuWAVELOOP(:,:,1);%.* exp(-1i*(TBkitobj.orbL*TBkitobj.Rm*(klist_tmp(end,:)-klist_tmp(1,:)).'));
        pmulist(:,ik_evolution) = TBkit.wancenter_1D(NuWAVELOOP);
    end
    % second Construct WANhamiltionian
    pmulist(pmulist<0) = pmulist(pmulist<0)+2*pi;
    nested_BFCAR(:,k) = mean(pmulist,2);
    pb.print(k,knum_evol,' ...');
end
pb.delete();
end