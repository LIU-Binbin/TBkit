function H_hr = Subsall(H_hr,mode,options)
% SUBSALL Complete substitution of all symbolic variables
%
%   H_HR = SUBSALL(H_HR,MODE) substitutes all symbolic variables
%
%   Inputs:
%       H_hr - HR object
%       mode - Substitution mode ('num','file','sym') [default: 'num']
%   Output:
%       H_hr - Modified HR object
%
%   Notes:
%       - 'num' mode converts to numeric
%       - 'file' mode loads from para.mat
%       - 'sym' mode performs symbolic simplification
%       - Checks for unresolved variables
arguments
    H_hr HR
    mode {mustBeMember(mode,{'num','file','add','sym'})} = 'num'
    options.OverlapTest = true;
end
if H_hr(1).overlap && options.OverlapTest
    H_hr(1) = Subsall(H_hr(1),mode,'OverlapTest',false);
    try
        H_hr(2) = Subsall(H_hr(2),mode,'OverlapTest',false);
    catch
    end
    return
end

if exist('para.mat','file') && strcmp(mode,'file')
    load('para.mat');
else
end
switch mode
    case {'num','file'}
        HcoeL_temp = subs(H_hr.HcoeL);
        symname = symvar(HcoeL_temp);
        if ~isempty(symname)
            for i = 1:length(symname)
                warning('this varible:');
            end
            disp(symname);
            error('please introduce the varible value above or consider use sym mode');
        else
            if isempty(HcoeL_temp)
            else
                H_hr.HnumL = double(HcoeL_temp);
                H_hr.HcoeL = sym([]);
            end
            H_hr.num = true;
            H_hr.coe = false;
        end
    case 'sym'
        H_hr.HcoeL = subs(H_hr.HcoeL);
    case 'add'
        H_hr.HnumL = H_hr.HnumL + double(subs(H_hr.HcoeL));
    otherwise
end
end
