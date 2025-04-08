function H_htrig = Subsall(H_htrig,mode,AddtionVar)
arguments
H_htrig
mode {mustBeMember(mode,{'num','file','gen','para','sym','slab','disk'})}= 'num';
AddtionVar = sym([]);
end
if nargin <2
mode = 'num';
end
if nargin < 3
AddtionVar = sym([]);
end
if exist('para.mat','file') && strcmp(mode,'file')
load('para.mat');
else
end
if ~isempty(AddtionVar)
varlist = [ H_htrig.VarsSeqLcart(1:H_htrig.Dim),AddtionVar];
else
varlist = H_htrig.VarsSeqLcart(1:H_htrig.Dim);
end
H_htrig.HcoeL = subs( H_htrig.HcoeL);
H_htrig.HsymL_trig = subs(H_htrig.HsymL_trig);
H_htrig.HsymL_coeL  = subs(H_htrig.HsymL_coeL);
H_htrig.HsymL_trig_bk = subs(H_htrig.HsymL_trig_bk);
H_htrig.Htrig_num  = subs(H_htrig.Htrig_sym);
H_htrig.Hsym = H_htrig.Htrig_num;
SymvarListInHcoeL = symvar(H_htrig.HcoeL);
switch mode
case {'num','file','gen','para'}
HcoeL_temp = H_htrig.HcoeL;
if ~isempty(SymvarListInHcoeL)
for i = 1:length(SymvarListInHcoeL)
fprintf('this varible: %s is the extra parameter\n',string(SymvarListInHcoeL(i)));
if isempty(find(varlist == SymvarListInHcoeL(i), 1))
try
HcoeL_temp = subs(HcoeL_temp,SymvarListInHcoeL(i),evalin('base',string(SymvarListInHcoeL(i))));
catch
warning('please subs this varible, Subsall must return numerical obj!\n');
end
end
end
else
H_htrig.HnumL = double(HcoeL_temp);
end
if strcmp(mode,'para')
H_htrig.HcoeL = subs(H_htrig.HcoeL);
return;
else
H_htrig.HnumL = double(HcoeL_temp);
H_htrig.HcoeL = sym([]);
end
H_htrig.num = true;
H_htrig.coe = false;
if strcmp(H_htrig.Type,'mat') ||strcmp(H_htrig.Type,'list')
H_htrig.HsymL_numL  = double(H_htrig.HsymL_coeL);
H_htrig.HsymL_coeL = sym([]);
return;
end
if strcmp(mode,'gen')
elseif ~strcmp(H_htrig.Type,'slab')
H_htrig.Htrig_num  = subs(H_htrig.Htrig_sym);
if H_htrig.Basis_num >500
H_htrig = sparse(H_htrig);
else
end
H_htrig.HcoeL =[];
else
if H_htrig.Basis_num >500
else
end
end
case {'sym','slab','disk'}
H_htrig.HcoeL = subs(H_htrig.HcoeL);
end
end
