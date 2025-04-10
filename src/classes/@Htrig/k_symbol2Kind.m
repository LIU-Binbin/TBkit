% K_SYMBOL2KIND Map a k_symbol to its corresponding kind in the Htrig object.
%
% SYNTAX:
%   Kind = k_symbol2Kind(H_htrig, k_symbol)
%
% DESCRIPTION:
%   This function finds the index (Kind) corresponding to the given k_symbol within the
%   Htrig object's symbolic hopping term list. The search method depends on the Type of 
%   the Htrig object:
%
%     - For 'exp' and 'sincos' types, the function converts the k_symbol and the elements 
%       of HsymL_trig to strings and then finds the matching index.
%
%     - For 'mat' and 'list' types, if k_symbol is symbolic, it uses ismember to find the 
%       row index in HsymL_coeL; otherwise, it searches in HsymL_numL.
%
%     - For any other types, it defaults to comparing string representations with HsymL_trig.
%
% INPUTS:
%   H_htrig  - A Htrig object.
%   k_symbol - A symbolic expression representing a k-space symbol.
%
% OUTPUT:
%   Kind     - The index (or indices) where k_symbol matches the corresponding entry in the 
%              Htrig object's hopping term list. If no match is found, Kind is empty.
%
% EXAMPLE:
%   % Given a Htrig object H and a symbolic k_symbol:
%   Kind = k_symbol2Kind(H, k_symbol);
%
function Kind = k_symbol2Kind(H_htrig, k_symbol)
    switch H_htrig.Type
        case {'exp','sincos'}
            str_2_compare = string(k_symbol);
            Kind = find(string(H_htrig.HsymL_trig) == str_2_compare);
        case {'mat','list'}
            if isa(k_symbol, 'sym')
                [~, Kind] = ismember(k_symbol, H_htrig.HsymL_coeL, 'rows');
            else
                [~, Kind] = ismember(k_symbol, H_htrig.HsymL_numL, 'rows');
            end
        otherwise
            str_2_compare = string(k_symbol);
            Kind = find(string(H_htrig.HsymL_trig) == str_2_compare);
    end
end

