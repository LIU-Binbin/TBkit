function H_hr = Subsall(H_hr,mode)
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
if nargin <2
mode = 'num';
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
if H_hr.overlap
ScoeL_temp = subs(H_hr.ScoeL);
symname = symvar(ScoeL_temp);
if ~isempty(symname)
for i = 1:length(symname)
warning('this varible:');
end
disp(symname);
error('please introduce the varible value above or consider use sym mode');
else
H_hr.SnumL = double(ScoeL_temp);
if ~strcmp(H_hr.Type,'list')
H_hr.num = true;
H_hr.coe = false;
end
end
end
case 'sym'
H_hr.HcoeL = subs(H_hr.HcoeL);
otherwise
end
end
