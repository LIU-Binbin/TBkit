function H_hk = Subsall(H_hk,mode)
%SUBSALL Complete symbolic substitution in HK object
%   H_HK = SUBSALL(H_HK, MODE) performs complete symbolic substitution.
%
%   Inputs:
%       H_hk - Input HK object
%       mode - Substitution mode:
%              'num'  - Convert to numeric (default)
%              'file' - Load substitutions from para.mat
%              'gen'  - Generic substitution
%              'sym'  - Symbolic substitution
%
%   Output:
%       H_hk - HK object after substitution
%
%   Note:
%       Requires all symbolic variables to be defined for numeric conversion
if nargin <2
    mode = 'num';
end
if exist('para.mat','file') && strcmp(mode,'file')
    load('para.mat');
else
end
switch mode
    case {'num','file','gen'}
        HcoeL_temp = subs(H_hk.HcoeL);
        symname = symvar(HcoeL_temp);
        if ~isempty(symname)
            for i = 1:length(symname)
                warning('this varible:
            end
            disp(symname);
            error('please introduce the varible value above');
        else
            H_hk.HnumL = double(HcoeL_temp);
        end
        H_hk.Hk_num  = subs(H_hk.Hk_sym);
        H_hk.Hsym  = H_hk.Hk_num ;
    case 'sym'
        H_hk.HcoeL = subs(H_hk.HcoeL);
        H_hk.Trig_to_save = subs(H_hk.Trig_to_save);
end
end
