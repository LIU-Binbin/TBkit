function sym_term = standardize_sym(sym_term)
%STANDARDIZE_SYM Normalize k-space variable notation
%
% Syntax:
%   std_term = standardize_sym(k_term)
%
% Input:
%   sym_term - Symbolic k-space expression
%
% Output:
%   sym_term - Standardized expression
%
% Description:
%   Ensures consistent variable naming by converting:
%   - All x-component variants (k_x, kX, Kx, etc.) → 'k_x'
%   - All y-component variants → 'k_y'
%   - All z-component variants → 'k_z'
%
% Note:
%   Essential for reliable term matching/comparison
%   Preserves mathematical structure while normalizing symbols
%
% Example:
%   std_term = standardize_sym('Kx*ky^2') % Returns 'k_x*k_y^2'
str_tmp = string(sym_term);
str_tmp = replace(str_tmp,{'k_x','k_X','K_X','K_x','kx','kX','KX','Kx'},'k_x');
str_tmp = replace(str_tmp,{'k_y','k_Y','K_y','K_Y','ky','kY','Ky','KY'},'k_y');
str_tmp = replace(str_tmp,{'k_z','k_Z','K_z','K_Z','kz','kZ','Kz','KZ'},'k_z');
sym_term = str2sym(str_tmp);
end
