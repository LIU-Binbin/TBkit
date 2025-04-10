function sym_coe_list = coeff_extract(sym_term)
% COEFF_EXTRACT Extracts variable components from symbolic product terms.
%   sym_coe_list = COEFF_EXTRACT(sym_term) decomposes the input symbolic 
%   product into subcomponents containing specified variables (x/y/z) and 
%   returns them as a symbolic array. Note: Coefficient factors are not 
%   preserved in current implementation.
%
%   Inputs:
%       sym_term - Symbolic expression in product form (e.g., 2*x*sin(y))
%
%   Outputs:
%       sym_coe_list - Symbolic array of components containing variables 
%                      from the product terms. Non-variable factors are 
%                      excluded.
%
%   Processing Steps:
%       1. Split input term into multiplicative components using '*' delimiter
%       2. Identify components containing x/y/z variables
%       3. Convert qualified components to symbolic form
%
%   Example:
%       term = sym('3*a*x*cos(y)');
%       components = coeff_extract(term) 
%       % Returns: [x, cos(y)] (a and constant factors are excluded)
%
%   Note: 
%       - Requires Htrig.strcontain helper function for variable detection
%       - Coefficients are currently not extracted (implementation limitation)
%       - Pre-allocation issues may exist in loop implementation
%
%   See also: sym, str2sym, Htrig.strcontain

% Split symbolic term into string components
str_list_tmp = strsplit(string(sym_term),'*');

% Initialize storage (Note: pre-allocation recommended for production code)
sym_coe = sym(1); % Placeholder coefficient (not utilized in current logic)
sym_coe_list = sym([]); % Output array initialization

% Component filtering loop
for i = 1:length(str_list_tmp)
    % Check for variable presence using helper function
    if Htrig.strcontain(str_list_tmp{i}, ['x','y','z'])
        % Store variable-containing components
        sym_coe_list(end+1) = sym_coe * str2sym(str_list_tmp{i});
    end
end
end