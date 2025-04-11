function H_htrig2 = subs(H_htrig, varargin)
%SUBS  Substitutes variables in the symbolic Hamiltonian and returns a modified Htrig object.
%
%   H_htrig2 = SUBS(H_htrig) substitutes all symbolic variables in the Hamiltonian terms 
%   of the input Htrig object (H_htrig) and returns the updated Htrig object (H_htrig2).
%   
%   H_htrig2 = SUBS(H_htrig, var1) substitutes the variable `var1` for all instances 
%   of the symbolic variables in the Hamiltonian terms of the input Htrig object.
%
%   H_htrig2 = SUBS(H_htrig, var1, var2) substitutes the variable `var1` for the first 
%   variable and `var2` for the second variable in the Hamiltonian terms of the input Htrig object.
%
%   Inputs:
%       H_htrig    - An instance of the Htrig class containing the symbolic Hamiltonian terms.
%       varargin   - A variable-length argument list containing the variables to substitute.
%                   If no arguments are provided, all symbolic variables in the Hamiltonian 
%                   will be substituted with their simplified expressions. If one or two variables 
%                   are provided, they will be used to perform the substitution in the Hamiltonian.
%
%   Outputs:
%       H_htrig2   - An updated Htrig object with substituted Hamiltonian terms.
%
%   Behavior:
%       - The function first checks the number of input arguments to determine which symbolic variables 
%         need to be substituted into the Hamiltonian.
%       - The Hamiltonian terms (stored in `HsymL_trig`) are expanded and simplified after substitution.
%       - The function iterates over the Hamiltonian kinds, splitting each symbolic term into its coefficients 
%         and associated symbolic variables.
%       - The resulting coefficients and variables are then organized into appropriate cell arrays.
%       - The Htrig object is updated with the new variables and coefficients, and returned as the output.
%
%   Example:
%       syms kx ky kz real
%       H_htrig = Htrig();                  % Create symbolic Htrig object
%       H_htrig2 = subs(H_htrig, kx, 1);    % Substitute kx = 1 into the Hamiltonian
%
%   See also: Htrig, split_sym_eq, simplify, subs, expand

    HsymL_trig_tmp = H_htrig.HsymL_trig;
    switch length(varargin)
        case 0
            HsymL_trig_tmp = expand(simplify(subs(HsymL_trig_tmp)));
        case 1
            HsymL_trig_tmp = expand(simplify(subs(HsymL_trig_tmp, varargin{1})));
        case 2
            HsymL_trig_tmp = expand(simplify(subs(HsymL_trig_tmp, varargin{1}, varargin{2})));
    end
    H_htrig2 = H_htrig;
    H_htrig2.HcoeL = sym([]);
    H_htrig2.HnumL = [];
    H_htrig2.HsymL_trig = sym([]);
    count = 0;
    for i = 1:H_htrig.Kinds
        [coeff_trig, symvar_list_trig, H_htrig2] = split_sym_eq(H_htrig2, HsymL_trig_tmp(i));
        for j = 1:numel(coeff_trig)
            count = count + 1;
            k_cell{count} = symvar_list_trig(j);
            mat_cell{count} = H_htrig.HcoeL(:,:,i);
            Var_cell{count} = coeff_trig(j);
        end
    end
    H_htrig2 = H_htrig2.setup(Var_cell, k_cell, mat_cell);
end
