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
%       - Handles both coefficients and overlaps
%       - Returns substitution equivalences
%       - Automatically simplifies results
if strcmp(varargin{end},'all')
all_mode = true;
nargin_check = nargin-1;
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
case 5
if strcmp(H_hr.Type , 'list')
H_hr.HcoeL = expand(simplify(subs(H_hr.HcoeL,varargin{1},varargin{2})));
H_hr.ScoeL = expand(simplify(subs(H_hr.ScoeL,varargin{3},varargin{4})));
else
H_hr.HcoeL = simplify(subs(H_hr.HcoeL,varargin{1},varargin{2}));
H_hr.ScoeL = simplify(subs(H_hr.ScoeL,varargin{3},varargin{4}));
end
EQL{1} =(varargin{1}==vpa(varargin{2}));
EQL{2} =(varargin{3}==vpa(varargin{3}));
end
if  isempty(H_hr.symvar_list) && ~H_hr.overlap
H_hr = H_hr.Subsall();
elseif isempty(H_hr.symvar_list) &&  H_hr.overlap && isempty(symvar(H_hr.ScoeL))
H_hr = H_hr.Subsall();
else
end
end
