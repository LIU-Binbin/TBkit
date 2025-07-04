function [BCCAR,Grid,BC_WAVECAR,klist_r_plot] = BC_2D(TBkitobj,optionsK,options,optionsPlot)
%BC_2D Calculate 2D Berry curvature distribution
%
%   Syntax:
%       [BCCAR,Grid,BC_WAVECAR,klist_r_plot] = BC_2D(TBkitobj,optionsK,options,optionsPlot)
%
%   Description:
%       Computes Berry curvature distribution over a 2D k-mesh for a
%       given tight-binding model. Supports both numerical and symbolic
%       Hamiltonians.
%
%   Inputs:
%       TBkitobj   - Tight-binding model object (HR, Htrig, or HK)
%       optionsK   - k-mesh options structure:
%                    knum1/2    - Grid dimensions
%                    kstart     - Starting k-point
%                    kdir1/2    - Grid directions
%                    cartesian  - Cartesian coordinates flag
%                    dir_seq    - Direction sequence
%                    dir_start  - Starting coordinate system
%       options    - Calculation options:
%                    BAND_index    - Band indices to include
%                    fig/plot      - Visualization flags
%                    Oper          - Operator for projection
%                    subband       - Subband indices
%                    sum           - Sum over bands flag
%                    ProjectionMethod - 'sign' or 'dynamicsign'
%       optionsPlot - Plotting options
%
%   Outputs:
%       BCCAR       - Berry curvature array
%       Grid        - k-mesh grid
%       BC_WAVECAR  - Wavefunctions used in calculation
%       klist_r_plot - k-points for plotting
arguments
    TBkitobj ;
    optionsK.knum1   = 51;
    optionsK.knum2   = 51;
    optionsK.kstart  = [-0.5,-0.5,0];
    optionsK.kdir1   = [1,0,0];
    optionsK.kdir2   = [0,1,0];
    optionsK.cartesian = false;
    optionsK.dir_seq = [1,2,3];
    optionsK.dir_start = 'kcar';
    options.BAND_index = [];
    options.fig = false;
    options.plot = false;
    options.Oper = [];
    options.subband = [];
    options.sum = true;
    options.ProjectionMethod {mustBeMember(options.ProjectionMethod,{'sign','dynamicsign'})}= 'sign';
    options.ProjectionStruct = struct('field','imag');
    optionsPlot.oneshot = true;
    optionsPlot.view = [0,90];
end
% --- nargin
switch class(TBkitobj)
    case {'Htrig','HK'}
        %
    case 'HR'
        TBkitobj.Basis_num = TBkitobj.WAN_NUM;
end
optionsKcell = namedargs2cell(optionsK);
[klist_cart,klist_frac,klist_r_plot,sizemesh,Gk_,Grid] = TBkit.kmesh2D(TBkitobj.Rm,optionsKcell{:},'full',true);
% individial chern number
if isempty(options.Oper)
    project = false;
    set_divide = 2;
    AOperU = [];
else
    project = true;
    set_divide = 4;
    OperObj = options.Oper; % only support one Oper now
    switch options.ProjectionMethod
        case 'sign'
            OperU = roundn(OperObj.U,-6);
            [AOperU,UOperU] = eig(OperU);
            if strcmp(options.ProjectionStruct.field,'imag')
                [AOperU,Ucheck] = sorteig(imag(UOperU),AOperU);
            elseif strcmp(options.ProjectionStruct.field,'real')
                [AOperU,Ucheck] = sorteig(real(UOperU),AOperU);
            end
        case 'dynamicsign'
            OperU = OperObj.U;
            switch class(TBkitobj)
                case {'Htrig','HK'}
                    klist = klist_cart;
                case 'HR'
                    klist = klist_frac;
            end
            AOperU = zeros(TBkitobj.Basis_num,TBkitobj.Basis_num,size(klist,1));
            for i = 1:size(klist,1)
                [AOperUtmp,UOperU] = eig(OperU(klist(i,1),klist(i,2),klist(i,3)));
                if strcmp(options.ProjectionStruct.field,'imag')
                    [AOperUtmp,Ucheck] = park.sorteig(imag(UOperU),AOperUtmp);
                elseif strcmp(options.ProjectionStruct.field,'real')
                    [AOperUtmp,Ucheck] = park.sorteig(real(UOperU),AOperUtmp);
                end
                AOperU(:,:,i) = AOperUtmp;
            end
        otherwise

    end
end
if isempty(options.BAND_index)
    switch class(TBkitobj)
        case {'Htrig','HK','HR'}
            BAND_index = 1:(TBkitobj.Basis_num/set_divide);
        case ''
            BAND_index = 1:(TBkitobj.Basis_num/set_divide);
    end
else
    BAND_index = options.BAND_index;
end
if project
    if isempty(options.subband)
        subband = 1:TBkitobj.Basis_num/2;
    else
        subband = options.subband;
    end
else
    subband = options.subband;
end
%
switch class(TBkitobj)
    case {'Htrig','HK'}
        [~,BC_WAVECAR] = TBkitobj.EIGENCAR_gen(...
            'klist',klist_cart);
    case 'HR'
        [~,BC_WAVECAR] = TBkitobj.EIGENCAR_gen(...
            'klist',klist_frac,...
            'convention','I','printmode',false, ...
            'Umat',AOperU,'subband',subband);
end
BC_WAVECAR = BC_WAVECAR(:,BAND_index,:);
%reshape(BC_WAVECAR(:,BAND_index,:),[TBkitobj.Basis_num options.knum1*options.knum2]);
if options.sum
    BCCAR = sum(TBkit.BerryCuvature_2D( BC_WAVECAR ,sizemesh),3);
else
    BCCAR = TBkit.BerryCuvature_2D( BC_WAVECAR ,sizemesh,'sum',false);
end
% remove the edge data

if options.plot
    knum1 = sizemesh(1);
    knum2 = sizemesh(2);
    dk_1 = (optionsK.kdir1)/knum1*Gk_;
    dk_2 = (optionsK.kdir2)/knum2*Gk_;
    optionsPlotcell = namedargs2cell(optionsPlot);
    if options.sum
        dSumL = BCCAR(:)/(2*pi);
    else
        dSumL = sum(BCCAR,3)/(2*pi);
        dSumL = dSumL(:);
    end
    TBkit.ShowSurfIntegral(TBkitobj.Gk,klist_r_plot,dk_1,dk_2,dSumL,optionsPlotcell{:});
end
end