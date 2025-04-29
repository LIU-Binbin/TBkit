function options_extra = FitOptionHelper(EIGENCAR_DFT,algorithm,options)
%FITOPTIONHELPER Configure options for band structure fitting
%
%   Syntax:
%       options_extra = FitOptionHelper(EIGENCAR_DFT,algorithm,options)
%
%   Description:
%       Sets up fitting parameters and weights for different band structure
%       comparison algorithms (insulator, metal, Dirac point search).
%
%   Inputs:
%       EIGENCAR_DFT - Reference DFT eigenvalues
%       algorithm    - Fitting algorithm type:
%                      'pure_comparison', 'dirac', 'insulator', 'metal'
%       options      - Configuration structure with fields:
%                      show, WAN_NUM, Noccu, klist_range, etc.
%
%   Output:
%       options_extra - Configured fitting options structure
%
%   See also: EIGENCAR_Value, Degeneracy_EIGENCAR
arguments
    EIGENCAR_DFT;
    algorithm {mustBeMember(algorithm,['pure_comparison','dirac','insulator','metal'])} = pure_comparison ;
    options.show logical = false;
    options.WAN_NUM {mustBeInteger} = -1;
    options.Noccu {mustBeInteger} = 1;
    options.klist_range  = ':';
    options.NBAND_range_DFT = 1:size(EIGENCAR_DFT,1)/2;
    options.NBAND_range = [];
    options.weight_list = [];
    options.KPOINTS = 'KPOINTS';
    options.highK = [1,size(EIGENCAR_DFT,2)];
    options.E_range = [];
    options.GapThreshold = 0.1;
    options.Node_punishment = 3;
    options.Node_insect_punishment = 0.5;
    options.Dirac_area = [];
    options.Dirac_width = 3;
    options.Dirac_punishment = 5;
    options.Gap_punishment = 5;
    options.Gapless_punishment = 5;
end
if options.WAN_NUM  == -1
    options_extra.WAN_NUM = options.Noccu*2;
else
    options_extra.WAN_NUM = options.WAN_NUM;
end
%
if isempty(options.NBAND_range )
    options_extra.NBAND_range = 1:options.Noccu;
else
    options_extra.NBAND_range = options.NBAND_range;
end
%
if isempty(options.E_range)
    options_extra.E_range  = [min(min(EIGENCAR_DFT(options_extra.NBAND_range,:))),...
        max(max(EIGENCAR_DFT(options_extra.NBAND_range,:)))];
else
    options_extra.E_range  = options.E_range;
end
%
%设定拟合范围：
options_extra.NBAND_range_DFT = options.NBAND_range_DFT;
options_extra.klist_range = options.klist_range;
%
options_extra.Noccu = options.Noccu;
options_extra.GapThreshold = options.GapThreshold;
options_extra.highK = options.highK;
%
options_extra.Node_punishment = options.Node_punishment;
options_extra.Node_insect_punishment = options.Node_insect_punishment;
options_extra.Dirac_punishment = options.Dirac_punishment;
options_extra.Gap_punishment = options.Gap_punishment;
options_extra.Gapless_punishment = options.Gapless_punishment;
%
if strcmp(algorithm,'dirac')
    if isempty(options.Dirac_area)
        Dirac_area = find(abs(EIGENCAR_DFT(options_extra.Noccu+1,:) - EIGENCAR_DFT(options_extra.Noccu,:))...
            <options_extra.GapThreshold);
        DiracA = [];
        for i = Dirac_area
            DiracA = [DiracA,i-options.Dirac_width : i+options.Dirac_width];
        end
        options_extra.Dirac_area = unique(DiracA);
    else
        options_extra.Dirac_area = options.Dirac_area;
    end
end
%默认大小和斜率等权
if isempty(options.weight_list)
    switch algorithm
        case 'pure_comparison'
            options_extra.weight_list = [1,1];
        case 'dirac'
            options_extra.weight_list = [1,1,1,1,1,1];
        case 'insulator'
            options_extra.weight_list = [1,1,1,1,1,1];
        case 'metal'
            options_extra.weight_list = [1,1,1,1,1,1];
        otherwise
            options_extra.weight_list = [1,1,1,1,1,1];
    end
else
    options_extra.weight_list = options.weight_list;
    % disp([';;']);
end
end