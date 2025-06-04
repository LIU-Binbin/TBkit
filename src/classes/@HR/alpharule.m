function H_hr = alpharule(H_hr, level_cut, mode, options)
%ALPHARULE Apply parameter scaling rules to tight-binding Hamiltonian
%   This function implements distance-dependent parameter scaling rules for
%   Slater-Koster parameters in tight-binding models. Supports multiple scaling
%   modes with automatic or manual reference distance (Rd) selection.
%
%   Inputs:
%       H_hr(1)       - HR object containing Hamiltonian data
%       level_cut  - Maximum neighbor shell to process (-1=auto detect)
%       mode       - Scaling rule mode:
%                    0: No operation
%                    1: Uniform alpha scaling
%                    2: Per-parameter alpha scaling
%       options.Rd - Reference distance(s) for scaling (default=-1 auto)
%       options.silence - Suppress information output (default=true)
%
%   Output:
%       H_hr(1)       - Modified HR object with parameter substitution rules
%
%   Features:
%       - Handles both Hamiltonian (HcoeL) and overlap (ScoeL) matrices
%       - Supports automatic reference distance detection from nearest neighbors
%       - Generates exponential scaling rules: V(R) = V(Rd)*exp(-(R-Rd)/Rd*alpha)
%       - Maintains symbolic/numeric consistency based on input parameters

% Argument validation with type checking
arguments
    H_hr HR;
    level_cut double {mustBeInteger} = -1;
    mode double = 0;
    options.Rd double = -1;
    options.silence logical = true;
end

% Early return for mode 0
if mode == 0
    return;
end

% Auto-detect neighbor shell count if needed
if level_cut == -1
    level_cut = length(H_hr(1).Rnn);
end

% Get neighbor distance information
Rnn = H_hr(1).nn_information(options.silence);

% Initialize base parameter templates
base_hamiltonian = ["VssS","VspS","VsdS","VppS","VpdS","VppP","VpdP","VddS","VddP","VddD"];
base_overlap = ["SssS","SspS","SsdS","SppS","SpdS","SppP","SpdP","SddS","SddP","SddD"];

% Find minimum-order parameters in symbolic expressions
[base_symvar, base_idx] = process_base_parameters(H_hr(1).symvar_list, base_hamiltonian, level_cut);

[base_symvar_S, base_idx_S] = process_base_parameters(symvar(H_hr(2).HcoeL), base_overlap, level_cut);


% Generate substitution rules for different shells
[V_subs, S_subs] = generate_substitution_rules(mode, level_cut, base_symvar, base_idx,base_symvar_S, base_idx_S, Rnn, options);

% Apply substitutions to Hamiltonian
EQ = [V_subs.names ]==[ V_subs.values];
disp(EQ);
H_hr(1).HcoeL = subs(H_hr(1).HcoeL, V_subs.names, V_subs.values);
if H_hr(1).overlap
    EQ2 = [S_subs.names ]==[ S_subs.values];
    disp(EQ2);
    H_hr(2).HcoeL = subs(H_hr(2).HcoeL, S_subs.names, S_subs.values);
end

% Nested helper functions
    function [base_sym, base_num] = process_base_parameters(sym_list, base_names, max_level)
        % Find minimum existing parameters for each base type
        base_sym = sym([]);
        base_num = [];
        for name = base_names
            for level = 1:max_level
                param = sym(sprintf("%s_%d", name, level));
                if ismember(string(param), string(sym_list))
                    base_sym(end+1) = param;
                    base_num(end+1) = level;
                    break;
                end
            end
        end
    end

    function [H_subs,S_subs] = generate_substitution_rules(mode, max_shell, base_sym, base_idx, base_symvar_S, base_idx_S,Rnn, opts)
        % Core substitution rule generator
        H_subs = struct('names',[], 'values',[]);
        S_subs = struct('names',[], 'values',[]);
        
        switch mode
            case 1 % Uniform alpha scaling
                alpha = sym('alpha');
                H_subs = create_uniform_subs(base_sym, base_idx, Rnn, max_shell, alpha, opts);
                if H_hr(1).overlap
                    alpha_S = sym('alpha_S');
                    S_subs = create_uniform_subs(base_symvar_S, base_idx_S, Rnn, max_shell, alpha_S, opts);
                end
                
            case 2 % Per-parameter alpha scaling
                H_subs = create_per_param_subs(base_sym, base_idx, Rnn, max_shell, opts);
                if H_hr(1).overlap
                    S_subs = create_per_param_subs(base_symvar_S, base_idx_S, Rnn, max_shell, opts);
                end
        end
    end

    function subs = create_uniform_subs(base_sym, base_idx, Rnn, max_shell, alpha, opts)
        % Create uniform alpha substitution rules
        subs.names = [];
        subs.values = [];
        for i = 1:length(base_sym)
            ref_level = base_idx(i);
            for shell = (ref_level+1):max_shell
                param_name = sym(sprintf("%s_%d", strrep(char(base_sym(i)), '_1', ''), shell));
                if opts.Rd == -1
                    Rd = Rnn(ref_level);
                else
                    Rd = opts.Rd(min(i, length(opts.Rd)));
                end
                exponent = -(Rnn(shell) - Rd)/Rd * alpha;
                subs.names(end+1) = param_name;
                subs.values(end+1) = base_sym(i) * exp(exponent);
            end
        end
    end

    function subs = create_per_param_subs(base_sym, base_idx, Rnn, max_shell, opts)
        % Create per-parameter alpha substitution rules
        subs.names = sym([]);
        subs.values = sym([]);
        for i = 1:length(base_sym)
            param_base = strsplit(char(base_sym(i)), '_');
            alpha_name = sym(sprintf("alpha_%s", strjoin(param_base(1:end-1), '_')));
            ref_level = base_idx(i);
            for shell = (ref_level+1):max_shell
                param_name = sym(sprintf("%s_%d", strjoin(param_base(1:end-1), '_'), shell));
                if opts.Rd == -1
                    Rd = Rnn(ref_level);
                else
                    Rd = opts.Rd(min(i, length(opts.Rd)));
                end
                exponent = -(Rnn(shell) - Rd)/Rd * alpha_name;
                subs.names(end+1) = param_name;
                subs.values(end+1) = base_sym(i) * exp(exponent);
            end
        end
    end
end