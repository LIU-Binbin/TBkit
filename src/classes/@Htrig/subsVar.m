function [H_htrig,EQL] = subsVar(H_htrig, varargin)
%SUBSVAR Substitutes values into the Hamiltonian and checks equality.
%
%   [H_htrig, EQL] = SUBSVAR(H_htrig, varargin) substitutes values into the Hamiltonian 
%   terms in Htrig (H_htrig) based on the given symbolic variables and their corresponding 
%   values in the input arguments. The function supports substitution in different modes, 
%   depending on the number of input arguments.
%
%   Inputs:
%       H_htrig    - An instance of the Htrig class, containing the Hamiltonian terms and 
%                    symbolic variables.
%       varargin   - Variable input arguments for substitution. These can include the symbolic 
%                    variables and their corresponding values to substitute.
%                    - If one input argument is provided, all symbolic variables are substituted.
%                    - If two input arguments are provided, the first is a list of symbolic variables, 
%                      and the second is their corresponding values.
%                    - If five input arguments are provided, substitutions are performed for both 
%                      the Hamiltonian (HcoeL) and additional terms (ScoeL).
%
%   Outputs:
%       H_htrig    - The updated Htrig object after performing the substitutions.
%       EQL        - A logical value or cell array containing the equality checks for the 
%                    symbolic variables after substitution. If all symbolic variables are substituted,
%                    it returns a logical indicating equality. For five inputs, it returns a cell array
%                    of logicals for each substitution.
%
%   Behavior:
%       - The function handles substitution in several modes:
%         1. With a single input argument, all symbolic variables are substituted.
%         2. With two input arguments, symbolic variables are substituted based on the input list and values.
%         3. With five input arguments, both the Hamiltonian and additional terms are substituted.
%       - The function also checks if symbolic variables are successfully substituted, returning the result 
%         in the variable EQL.
%
%   Example:
%       % Substitute all symbolic variables in the Hamiltonian
%       [H_htrig, EQL] = subsVar(H_htrig, value);
%
%       % Substitute specific symbolic variables in the Hamiltonian
%       [H_htrig, EQL] = subsVar(H_htrig, symvar_list, values);
%
%       % Substitute values in both the Hamiltonian and ScoeL
%       [H_htrig, EQL] = subsVar(H_htrig, symvar_list1, values1, symvar_list2, values2);
%
%   See also: Htrig, Subsall

    if strcmp(varargin{end},'all')
        all_mode = true;
        nargin_check = nargin-1;
    else
        all_mode = false;
        nargin_check = nargin;
    end
    
    switch nargin_check
        case 1
            H_htrig.HcoeL = expand(simplify(subs(H_htrig.HcoeL)));
        case 2
            SymVarL = H_htrig.symvar_list;
            if all_mode
                H_htrig.HcoeL = subs(H_htrig.HcoeL,SymVarL,varargin{1});
            else
                H_htrig.HcoeL = (simplify(subs(H_htrig.HcoeL,SymVarL,varargin{1})));
            end
            EQL = (SymVarL == vpa(varargin{1}));
        case 3
            if all_mode
                H_htrig.HcoeL = subs(H_htrig.HcoeL,varargin{1},varargin{2});
            else
                H_htrig.HcoeL = simplify(subs(H_htrig.HcoeL,varargin{1},varargin{2}));
            end
            EQL = (varargin{1} == vpa(varargin{2}));
        case 5
            H_htrig.HcoeL = simplify(subs(H_htrig.HcoeL,varargin{1},varargin{2}));
            H_htrig.ScoeL = simplify(subs(H_htrig.ScoeL,varargin{3},varargin{4}));
            EQL{1} = (varargin{1} == vpa(varargin{2}));
            EQL{2} = (varargin{3} == vpa(varargin{3}));
    end

    if isempty(H_htrig.symvar_list)
        H_htrig = H_htrig.Subsall();
    else
    end
end

