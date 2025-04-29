function [vector_list,Coeffs_list] = HstrL_classify(H_strL_tmp,R_struct,mode)
%HSTRL_CLASSIFY Classify Hamiltonian terms for real-space mapping
%   [v, c] = HstrL_classify(H_strL_tmp, R_struct, mode) prepares symbolic
%   Hamiltonian terms for real-space vector mapping.
%
%   Inputs:
%       H_strL_tmp - String or symbolic expression
%       R_struct   - Lattice parameter structure
%       mode       - Lattice type ('T','H','H2')
%
%   Outputs:
%       vector_list  - Hopping vectors
%       Coeffs_list  - Corresponding coefficients
%
%   Example:
%       [v, c] = HstrL_classify('x+y', struct('a',1), 'H');
if nargin <3
mode = 'T';
end
if strcmp(mode,'T')
symvar_list = symvar(str2sym(H_strL_tmp));
if isempty(symvar_list)
vector_list = [0,0,0];
Coeffs_list = 1;
elseif length(symvar_list) == 1
[vector_list,Coeffs_list] = HK.Hstr_mapping(H_strL_tmp,R_struct,mode);
else
[vector_list,Coeffs_list] = HK.Hstr_mapping(string(symvar_list(1)),R_struct,mode);
for i = 2:length(symvar_list)
[VL2,CL2] = HK.Hstr_mapping(string(symvar_list(i)),R_struct,mode);
[vector_list,Coeffs_list] = HK.VLCL_ltimes(vector_list,Coeffs_list,VL2,CL2);
end
end
elseif strcmp(mode,'H')
[vector_list,Coeffs_list] = HK.Hstr_mapping(H_strL_tmp,R_struct,mode);
elseif strcmp(mode,'H2')
[vector_list,Coeffs_list] = HK.Hstr_mapping(H_strL_tmp,R_struct,mode);
end
end
