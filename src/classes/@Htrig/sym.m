function Htrig_sym = sym(H_htrig, options)
%SYM Generates a symbolic expression of the Hamiltonian in a specified form.
%
%   Htrig_sym = SYM(H_htrig, options) generates the symbolic expression for the Hamiltonian 
%   stored in the Htrig object (H_htrig) and returns it in the specified format. The function
%   allows customization of the output type (e.g., 'exp', 'sincos', etc.) and can simplify the 
%   result based on the options provided.
%
%   Inputs:
%       H_htrig   - An instance of the Htrig class containing the Hamiltonian and related properties.
%       options   - A structure with optional fields:
%           .simple (logical) - If true, the output will be simplified further (default: false).
%           .Type (string) - Specifies the type of symbolic expression. Options are:
%                            'exp', 'sincos', 'sin', 'cos', 'cot', 'tan', 'sinh', 'cosh', 
%                            'tanh', or an empty string (default: based on H_htrig.Type).
%                            
%   Outputs:
%       Htrig_sym - The symbolic representation of the Hamiltonian in the specified type.
%
%   Behavior:
%       - The function chooses the default symbolic type based on the Htrig object’s type if 
%         the 'Type' option is not provided.
%       - If the 'simple' option is true, it further simplifies the output expression by 
%         dividing by a specific term from the Htrig object.
%       - The result is rewritten in the specified form (such as 'exp' or 'sincos') and simplified 
%         where applicable.
%
%   Example:
%       % Generate the symbolic expression for Htrig in 'exp' format
%       Htrig_sym = sym(H_htrig, struct('Type', 'exp'));
%
%       % Generate the simplified symbolic expression for Htrig in 'sincos' format
%       Htrig_sym = sym(H_htrig, struct('simple', true, 'Type', 'sincos'));
%
%   See also: Htrig, timtj_gen

    Htrig_sym = H_htrig.Htrig_sym;
    
    % Determine the type of symbolic expression to use
    if strcmp(options.Type,'')
        switch H_htrig.Type
            case {'sincos'}
                type = 'sincos';
            case {'exp','list','mat'}
                type = 'exp';
            otherwise
                type = 'sincos';
        end
    else
        type = options.Type;
    end
    
    % Simplify the expression if the 'simple' option is enabled
    if options.simple
        H_htrig = H_htrig.timtj_gen('sym');
        Htrig_sym = simplify(Htrig_sym ./ H_htrig.timtj{3}.');
    end
    
    % Rewrite and simplify the symbolic expression as per the chosen type
    Htrig_sym = simplify(rewrite(Htrig_sym, type), 'IgnoreAnalyticConstraints', true);
end

