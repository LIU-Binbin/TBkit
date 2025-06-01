function [H_hr,EQL] = subs(H_hr,varargin)
% SUBS Symbolic substitution for HR object
%
%   [H_HR,EQL] = SUBS(H_HR,VARARGIN) performs symbolic substitution
%
%   Inputs:
%       H_hr - HR object
%       varargin - Substitution arguments (see MATLAB subs)
%   Outputs:
%       H_hr - Modified HR object
%       EQL - Equality list
%
%   Notes:
%       - Supports 'all' mode for full substitution
%       - Returns substitution equivalences
%       - Automatically simplifies results

if strcmp(varargin{end},'all')
    all_mode = true;
    nargin_check = nargin-1;
elseif strcmp(varargin{end-1},'SubsIndexL')
    switch nargin-2
        case 2
            SubsIndexL = varargin{end};
            x = varargin{1};
            symvar_list = H_hr.symvar_list;
            H_hr.HnumL(SubsIndexL) = double(subs(H_hr.HcoeL,symvar_list,x));
            EQL = vpa(symvar_list == x);
        case 3
            SubsIndexL = varargin{end};
            x = varargin{2};
            symvar_list = varargin{1};
            %SubsHnumL = double(subs(FITobj.HcoeL,Varlist,parameters));
            H_hr.HnumL(SubsIndexL) = H_hr.HnumL(SubsIndexL) + double(subs(H_hr.HcoeL,symvar_list,x));
            EQL = vpa(symvar_list == x);
    end
    return;
else
    all_mode = false;
    nargin_check = nargin;
end
switch nargin_check
    case 1
        H_hr.HcoeL = expand(simplify(subs(H_hr.HcoeL)));
    case 2
        SymVarL = H_hr.symvar_list;
        if all_mode
            H_hr.HcoeL = subs(H_hr.HcoeL,SymVarL,varargin{1});
        else
            if strcmp(H_hr.Type , 'list')
                H_hr.HcoeL = expand(simplify(subs(H_hr.HcoeL,SymVarL,varargin{1})));
            else
                H_hr.HcoeL = (simplify(subs(H_hr.HcoeL,SymVarL,varargin{1})));
            end
        end
        EQL =(SymVarL==vpa(varargin{1}));
    case 3
        if all_mode
            H_hr.HcoeL = subs(H_hr.HcoeL,varargin{1},varargin{2});
        else
            if strcmp(H_hr.Type , 'list')
                H_hr.HcoeL = expand(simplify(subs(H_hr.HcoeL,varargin{1},varargin{2})));
            else
                H_hr.HcoeL = simplify(subs(H_hr.HcoeL,varargin{1},varargin{2}));
            end
        end
        EQL =(varargin{1}==vpa(varargin{2}));
end
if  isempty(H_hr.symvar_list)
    H_hr = H_hr.Subsall();
else
    
end
end
