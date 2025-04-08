function [H_htrig,EQL] = subsVar(H_htrig,varargin)
if strcmp(varargin{end},'all')
all_mode = true;
nargin_check = nargin-1;
else
all_mode = false;
nargin_check = nargin;
end
switch nargin_check
case 1
H_htrig.HcoeL = expand(simplify(subs(H_htrig.HcoeL)));
case 2
SymVarL = H_htrig.symvar_list;
if all_mode
H_htrig.HcoeL = subs(H_htrig.HcoeL,SymVarL,varargin{1});
else
H_htrig.HcoeL = (simplify(subs(H_htrig.HcoeL,SymVarL,varargin{1})));
end
EQL =(SymVarL==vpa(varargin{1}));
case 3
if all_mode
H_htrig.HcoeL = subs(H_htrig.HcoeL,varargin{1},varargin{2});
else
H_htrig.HcoeL = simplify(subs(H_htrig.HcoeL,varargin{1},varargin{2}));
end
EQL =(varargin{1}==vpa(varargin{2}));
case 5
H_htrig.HcoeL = simplify(subs(H_htrig.HcoeL,varargin{1},varargin{2}));
H_htrig.ScoeL = simplify(subs(H_htrig.ScoeL,varargin{3},varargin{4}));
EQL{1} =(varargin{1}==vpa(varargin{2}));
EQL{2} =(varargin{3}==vpa(varargin{3}));
end
if  isempty(H_htrig.symvar_list)
H_htrig = H_htrig.Subsall();
else
end
end
