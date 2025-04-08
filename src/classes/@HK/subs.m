function H_hk  = subs(H_hk,varargin)
switch length(varargin)
case 1
H_hk.HcoeL = subs(H_hk.HcoeL,varargin{1});
H_hk.Trig_to_save = subs(H_hk.Trig_to_save,varargin{1});
case 2
H_hk.HcoeL = subs(H_hk.HcoeL,varargin{1},varargin{2});
H_hk.Trig_to_save = subs(H_hk.Trig_to_save,varargin{1},varargin{2});
case 3
H_hk.HcoeL = subs(H_hk.HcoeL,varargin{1},varargin{2},varargin{3});
H_hk.Trig_to_save = subs(H_hk.Trig_to_save,varargin{1},varargin{2},varargin{3});
end
end
