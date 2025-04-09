function H_hk = Subsall(H_hk, mode)
%SUBSALL Complete parameter substitution
%
% Syntax:
%   Hk = Subsall(Hk)
%   Hk = Subsall(Hk, mode)
%
% Inputs:
%   Hk - Original Hamiltonian
%   mode - Substitution mode:
%       'num' - Numeric substitution (default)
%       'file' - Load values from para.mat
%       'sym' - Symbolic only
%       'gen' - Generic numeric
%
% Output:
%   Hk - Hamiltonian with substitutions
%
% Description:
%   Comprehensive substitution that:
%   1. Automatically detects remaining variables
%   2. Handles different substitution sources
%   3. Converts to numeric form when complete
%
% Note:
%   Throws error if unresolved variables remain
%   'file' mode requires para.mat with variables
%
% Example:
%   Hk_num = Hk.Subsall('num'); % Full numeric conversion
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
                warning('this varible: %s is not defind',string(symname(i)));
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
