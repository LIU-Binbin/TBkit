function Kind = k_symbol2Kind(H_hk,k_symbol)
%K_SYMBOL2KIND Map k-space symbol to term index
%
% Syntax:
%   Kind = k_symbol2Kind(H_hk,k_symbol)
%
% Inputs:
%   H_hk - HK object containing term definitions
%   k_symbol - String/sym representing k-space term (e.g., 'k_x^2')
%
% Output:
%   Kind - Index of matching term in HsymL
%
% Description:
%   Normalizes k-space expressions and matches against stored terms by:
%   1. Standardizing variable names (k_x, k_y, k_z)
%   2. Extracting coefficients
%   3. Comparing normalized forms with HsymL terms
%
% Note:
%   Returns empty if no match found
%   Case-insensitive for variable names
%
% Example:
%   idx = k_symbol2Kind(Hk, 'K_y^2') % Finds Y^2 term index
k_symbol = string(k_symbol);
symvar_list = symvar(str2sym(k_symbol));
for i = 1:length(symvar_list)
    str_tmp = string(symvar_list(i));
    switch str_tmp
        case {'k_x','k_X','K_X','K_x','kx','kX','KX','Kx','x','X'}
            k_symbol = strrep(k_symbol,str_tmp,'k_x');
        case {'k_y','k_Y','K_y','K_Y','ky','kY','Ky','KY','y','Y'}
            k_symbol = strrep(k_symbol,str_tmp,'k_y');
        case {'k_z','k_Z','K_z','K_Z','kz','kZ','Kz','KZ','z','Z'}
            k_symbol = strrep(k_symbol,str_tmp,'k_z');
        case {'k_w','k_W','K_w','K_W','kw','kW','Kw','KW'}
            k_symbol_str = strrep(k_symbol_str,str_tmp,'k_w');
    end
end
coeff_tmp = coeffs(str2sym(k_symbol));
str_2_compare = string((str2sym(k_symbol)/coeff_tmp));
Kind = find(string(H_hk.HsymL) == str_2_compare);
end
