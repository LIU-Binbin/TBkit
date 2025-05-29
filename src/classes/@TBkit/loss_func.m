function Value = loss_func(parameters,extra_parm,options)
%LOSS_FUNC Compute loss function for tight-binding parameter optimization
%
%   Syntax:
%       Value = loss_func(parameters,extra_parm,options)
%
%   Description:
%       Evaluates quality metric between DFT and model band structures
%       for parameter optimization. Supports multiple fitting methods.
%
%   Inputs:
%       parameters - Parameters to evaluate (table or double array)
%       extra_parm - Additional weighting parameter (default=0)
%       options    - Configuration structure with fields:
%                    mode, algorithm, DFTBAND, FITobj, etc.
%
%   Output:
%       Value - Computed loss metric
%
%   See also: EIGENCAR_Value, FitOptionHelper
arguments
    parameters;
    extra_parm double =0;
    options.mode = 'extra';
    options.algorithm  = 'pure_comparison';
    options.DFTBAND = 'EIGENCAR_DFT';
    options.FITobj = 'H_TBkit';
    options.Varlist = 'Varlist';
    options.Varlock = 'Varlock';
    options.namespace = 'base';
    options.extra = 'options_extra';
    options.show = false;
    options.fig = handle([]);
    options.ax = handle([]);
    options.KPOINTS = 'KPOINTS';
    options.SubsIndexL = []; 
end
%%%%%%%%%%%%%%
EIGENCAR_DFT = evalin(options.namespace,options.DFTBAND);
FITobj = evalin(options.namespace,options.FITobj);
if isempty(options.SubsIndexL)
    SubsIndexL = 1:numel(FITobj.HnumL);
else
    SubsIndexL = options.SubsIndexL;
end


try
    Varlist = evalin(options.namespace,options.Varlist);
catch
    try
        Varlist = FITobj.symvar_list;
    catch
        Varlist = [];
    end
end

%%%%%%%%%%%%%%%
if isa(parameters,'table')
    fitmethod = 'Bayes';
elseif isa(parameters,'double')
    if length(parameters) > 1
        fitmethod = 'NM';
    else
        fitmethod = 'Single';
    end
end
%%%%%%%%%%%%%%%
switch fitmethod
    case {'Single','NM'}
        if isa(FITobj,'function_handle')
            FITobj_n = FITobj(parameters);
            FITobj_n = FITobj_n <options.KPOINTS;
        else
            try
                SubsHnumL = double(subs(FITobj.HcoeL,Varlist,parameters));
            catch
                Varlist = FITobj.symvar_list;
                SubsHnumL = double(subs(FITobj.HcoeL,Varlist,parameters));
            end
            FITobj_n = FITobj;
            FITobj_n.HnumL(SubsIndexL) = SubsHnumL;
        end
        EIGENCAR_TBkit = FITobj_n.EIGENCAR_gen('printmode',false);
    case 'Bayes'
        %%%%%%%%%%%%%%%%%
        try
            Varlock = evalin(options.namespace,options.Varlock);
        catch

        end
        SubsHnumL = FITobj.HcoeL;
        %%%%%%%%%%%%%%%%%
        for i  = 1:length(Varlist)
            % Generate Field Names from Variables  dynamic fieldnames, or sometimes dynamic field names.
            %FITobj = FITobj.subs(Varlist(i),parameters.(string(Varlist(i))));
            try
                SubsHnumL = (subs(SubsHnumL,Varlist,parameters.(string(Varlist(i))) ));
            catch
                Varlist = FITobj.symvar_list;
                SubsHnumL = (subs(SubsHnumL,Varlist,parameters.(string(Varlist(i))) ));
            end
        end

        FITobj_n = FITobj;
        FITobj_n.HnumL(SubsIndexL) = SubsHnumL;
        EIGENCAR_TBkit = FITobj_n.EIGENCAR_gen('printmode',false);
    otherwise
        error('not be implemented');
end
Value = TBkit.EIGENCAR_Value(EIGENCAR_DFT,EIGENCAR_TBkit,extra_parm,...
    'mode',options.mode,...
    'algorithm',options.algorithm,...
    'namespace',options.namespace,...
    'extra',options.extra,...
    'show',options.show,...
    'fig',options.fig,...
    'ax',options.ax ...
    );
end