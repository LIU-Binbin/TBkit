function ValueTotal = EIGENCAR_Value(EIGENCAR_DFT,EIGENCAR,extra_parm,options)
%EIGENCAR_VALUE Compare eigenvalue arrays with multiple metrics
%
%   Syntax:
%       ValueTotal = EIGENCAR_Value(EIGENCAR_DFT,EIGENCAR,extra_parm,options)
%
%   Description:
%       Computes similarity metrics between DFT and model eigenvalue arrays
%       with configurable weighting schemes for different physical properties.
%
%   Inputs:
%       EIGENCAR_DFT - Reference DFT eigenvalues
%       EIGENCAR     - Model eigenvalues to evaluate
%       extra_parm   - Additional parameter for weighting
%       options      - Configuration structure with fields:
%                      mode, namespace, extra, algorithm, show, fig, ax
%
%   Output:
%       ValueTotal - Composite quality metric value
%
%   See also: FitOptionHelper, EIGENCAR2IMG
arguments
    EIGENCAR_DFT double;
    EIGENCAR double;
    extra_parm double = 0;
    options.mode = 'extra';
    options.namespace = 'base';
    options.extra = 'options_extra'
    options.algorithm  = 'pure_comparison';
    options.show = false;
    options.fig = handle([]);
    options.ax = handle([]);
end
%%%%%%%%%%%%
if  ~strcmp(options.mode,'extra')
    DATA1 = EIGENCAR_DFT;
    DATA2 = EIGENCAR;
    weight_list = ones(2,1);
    Echoose = ':';
elseif strcmp(options.mode,'extra')
    options_extra = evalin(options.namespace,options.extra);
    NBAND_range_DFT = options_extra.NBAND_range_DFT;
    NBAND_range = options_extra.NBAND_range;
    klist_range = options_extra.klist_range;
    try
        weight_list = options_extra.weight_list;
    catch
        weight_list = ones(2,1);
    end
    DATA1 = EIGENCAR_DFT(NBAND_range_DFT,klist_range);
    DATA2 = EIGENCAR(NBAND_range,klist_range);
    if isfield(options_extra, 'E_range')
        if ~isempty(options_extra.E_range)
            E_range =  options_extra.E_range;
        else
            E_range = [min(DATA1,[],'all'),max(DATA1,[],'all')];
        end
    else
        E_range = [min(DATA1,[],'all'),max(DATA1,[],'all')];
    end
    Echoose =DATA1 >= E_range(1)  & DATA1 <= E_range(2);
end
%%%%%%%%%%%%%%%%%%%%
% Alldata
count = 0;
% Direct minus
count = count +1 ;

Value(count) = sqrt(mean(mean(abs(DATA1(Echoose)-DATA2(Echoose))).^2));
% Direct diff minus
count = count +1 ;
Value(count) = sqrt(mean(mean(abs(diff(DATA1(Echoose).')-diff(DATA2(Echoose).')))));
if strcmp(options.algorithm,'pure_comparison')
    ValueTotal = sum(Value.*weight_list);
    return;
else
    Noccu = options_extra.Noccu;
    GapThreshold = options_extra.GapThreshold;
    E_range =  options_extra.E_range;
    highK = options_extra.highK;
    GapL = (abs(EIGENCAR(Noccu +1,:)-EIGENCAR(Noccu ,:)));...
        minGap = abs((min(EIGENCAR_DFT(Noccu +1,:)-EIGENCAR_DFT(Noccu ,:)) - min(GapL)));
    Value_Gap_Label =min(GapL)> GapThreshold;
    % Node Structure
    Node_punishment = options_extra.Node_punishment;
    Node_insect_punishment = options_extra.Node_insect_punishment;
    [DEGENCAR1,NODEINSECT1] = TBkit.Degeneracy_EIGENCAR(EIGENCAR_DFT(1:Noccu+1,:),highK,GapThreshold/2);
    [DEGENCAR2,NODEINSECT2] = TBkit.Degeneracy_EIGENCAR(EIGENCAR(1:Noccu+1,:),highK,GapThreshold/2);
    Value_Node = Node_punishment*sum(sum(abs(DEGENCAR1-DEGENCAR2))) ...
        + Node_insect_punishment*sum(sum(abs(NODEINSECT1-NODEINSECT2)));
    % Flash
    IM1 = TBkit.EIGENCAR2IMG(EIGENCAR_DFT,GapThreshold/2,E_range);
    IM2 = TBkit.EIGENCAR2IMG(EIGENCAR,GapThreshold/2,E_range);
    Value_ssimval = 1/ssim(IM1,IM2)-1;
    % High symmetry line energy
    E_highK_DFT = EIGENCAR_DFT(NBAND_range_DFT,highK);
    E_highK_TB = EIGENCAR(NBAND_range_DFT,highK);
    Value_E_HighK = mean(mean(abs(E_highK_DFT-E_highK_TB)));
end
if strcmp(options.algorithm,'dirac')
    Dirac_area = options_extra.Dirac_area;
    Dirac_punishment = options_extra.Dirac_punishment;
    Gapless_points = GapL<0.3;
    Value_Dirac = 0;
    for i = 1:length(Gapless_points )
        if Gapless_points(i) == 1 && ~ismember(i,Dirac_area)
            Value_Dirac = Value_Dirac+ Dirac_punishment;
        end
    end
end
switch options.algorithm
    case 'dirac'
        count = count +1 ;
        Value(count) = Value_E_HighK;
        ValueTotal = (1-extra_parm)*(Value_ssimval^2)*sum(Value.*weight_list(1:3));
        %
        ValueTotal = ValueTotal...
            +extra_parm * Value_Gap_Label * weight_list(4)...
            +extra_parm *(Value_ssimval^3)  *(Value_Node+1) * weight_list(5)...
            +extra_parm * Value_Dirac * weight_list(6);
    case 'insulator'
        Gapless_punishment = options_extra.Gapless_punishment;
        count = count +1 ;
        Value(count) = Value_E_HighK;
        ValueTotal = (1-extra_parm) * (Value_ssimval^2) * sum(Value.*weight_list(1:3));
        %
        ValueTotal = ValueTotal...
            +extra_parm * minGap * weight_list(4)...
            +extra_parm * (Value_ssimval^3) * (Value_Node+1) * weight_list(5)...
            +extra_parm * Gapless_punishment * ~Value_Gap_Label * weight_list(6);
    case 'metal'
        Gap_punishment = options_extra.Gap_punishment;
        count = count +1 ;
        Value(count) = Value_E_HighK;
        ValueTotal = (1-extra_parm)*(Value_ssimval^2)*sum(Value.*weight_list(1:3));
        %
        ValueTotal = ValueTotal...
            +extra_parm * maxGap * weight_list(4)...
            +extra_parm * (Value_ssimval^3) * (Value_Node+1) * weight_list(5)...
            +extra_parm * Gap_punishment * Value_Gap_Label * weight_list(6);
    otherwise
end

%ValueTotal = ValueTotal^2;

end