function [vector_list,Coeffs_list] = HstrL_classify(H_strL_tmp,R_struct,mode)
%HSTRL_CLASSIFY Classify and map Hamiltonian string terms
%
% Syntax:
%   [v,c] = HstrL_classify(H_strL_tmp,R_struct)
%   [v,c] = HstrL_classify(H_strL_tmp,R_struct,mode)
%
% Inputs:
%   H_strL_tmp - String or symbolic term to classify
%   R_struct - Lattice constant structure
%   mode - Coordinate system ('T','H','H2')
%
% Outputs:
%   vector_list - Hopping vectors
%   Coeffs_list - Corresponding coefficients
%
% Description:
%   Wrapper for Hstr_mapping that handles:
%   - Constant terms (no variables)
%   - Single variable terms
%   - Multi-variable terms through recursive decomposition
%   - Different coordinate system modes
%
% Algorithm:
%   1. Extracts symbolic variables from term
%   2. For multi-variable terms, recursively combines mappings using VLCL_ltimes
%   3. Delegates actual mapping to Hstr_mapping
%
% See also:
%   Hstr_mapping, VLCL_ltimes
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
