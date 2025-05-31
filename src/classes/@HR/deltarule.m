function H_hr = deltarule(H_hr, level_cut, mode, options)
% DELTARULE Applies exponential decay substitution to Hamiltonian and overlap matrix coefficients
%   This function modifies the symbolic coefficients of an HR object's Hamiltonian and overlap matrix
%   using a delta rule substitution. It replaces higher-order interaction terms with exponentially
%   decaying functions of lower-order terms, based on nearest-neighbor distances.
%
% INPUTS:
%   H_hr       - HR object containing Hamiltonian (HcoeL) and overlap (ScoeL) symbolic expressions
%   level_cut  - Maximum interaction order to consider (default: 1, use -1 for full interaction range)
%   mode       - Substitution mode (0: no substitution, 1: single delta parameter, 2: per-term delta)
%   options.Rd - Custom decay reference distance(s) (default: -1 uses H_hr.Rnn)
%
% OUTPUTS:
%   H_hr       - Modified HR object with updated symbolic coefficients
%
% NOTES:
%   - mode 1 uses a single delta parameter for all substitutions
%   - mode 2 uses unique delta parameters for each interaction type
%   - level_cut determines the maximum interaction order (e.g., 3 includes 1st, 2nd, 3rd neighbors)
%   - Valid interaction types include VssS, VspS, VppS, etc., for Hamiltonian
%   - SssS, SspS, SppS, etc., are supported for overlap matrix when H_hr.overlap is true
%   - The substitution formula is: V_n = V_1 * exp(-(Rnn(j) - Rd)/delta)

arguments
    H_hr HR;
    level_cut double{mustBeInteger, mustBeNonnegative} = 1;
    mode double{mustBeMember(mode, [0, 1, 2])} = 0;
    options.Rd = -1;
end

% Handle special cases
if mode == 0
    return;
end

if level_cut == -1
    level_cut = length(H_hr.Rnn);
end

% Extract nearest-neighbor distances
Rnn = H_hr(1).nn_information();

% Define interaction types and symbolic variable patterns
base_string = ["VssS", "VspS", "VsdS", "VppS", "VpdS", "VppP", "VpdP", "VddS", "VddP", "VddD"];
strvar_list = string(H_hr(1).symvar_list);

% Initialize substitution parameters for Hamiltonian
base_symvar = sym([]);
base_num = [];
count = 0;

% Identify lowest-order terms for each interaction type
for ibase = base_string
    for i = 1:level_cut
        term = string(ibase) + "_" + num2str(i);
        if any(contains(strvar_list, term))
            count = count + 1;
            base_symvar(count) = sym(term, 'real');
            base_num(count) = i;
            break;
        end
    end
end

% Handle overlap matrix if present
if H_hr(1).overlap
    base_string_S = ["SssS", "SspS", "SsdS", "SppS", "SpdS", "SppP", "SpdP", "SddS", "SddP", "SddD"];
    symvar_list_S = string(symvar(H_hr(2).HcoeL));

    base_symvar_S = sym([]);
    base_num_S = [];
    count_S = 0;

    for ibase_S = base_string_S
        for i = 1:level_cut
            term = string(ibase_S) + "_" + num2str(i);
            if any(contains(symvar_list_S, term))
                count_S = count_S + 1;
                base_symvar_S(count_S) = sym(term, 'real');
                base_num_S(count_S) = i;
                break;
            end
        end
    end
end

% Configure decay reference distances
if options.Rd == -1
    RdL = Rnn;
else
    if length(options.Rd) > 1
        RdL = repmat(options.Rd, [1, count]);
    else
        RdL = options.Rd;
    end
    base_num = 1:count; % Reset to use all identified terms
    if H_hr(1).overlap
        base_num_S = 1:count_S;
    end
end

% Apply substitution rules based on mode
switch mode
    case 1
        % Single delta parameter for all interactions
        delta = sym('delta', 'real');

        % Substitute Hamiltonian terms
        for j = 2:level_cut
            for i = 1:count
                V_n = string(base_symvar(i)) + "_" + num2str(j);
                if any(contains(strvar_list, V_n))
                    coeff = (Rnn(j) - RdL(base_num(i)));
                    subs_expr = base_symvar(i) * exp(-coeff / delta);
                    H_hr(1).HcoeL = subs(H_hr(1).HcoeL, sym(V_n), subs_expr);
                end
            end
        end

        % Substitute overlap terms if present
        if H_hr(1).overlap
            delta_S = sym('delta__2', 'real');
            for j = 2:level_cut
                for i = 1:count_S
                    S_n = string(base_symvar_S(i)) + "_" + num2str(j);
                    if any(contains(symvar_list_S, S_n))
                        coeff = (Rnn(j) - RdL(base_num_S(i)));
                        subs_expr = base_symvar_S(i) * exp(-coeff / delta_S);
                        H_hr(2).HcoeL = subs(H_hr(2).HcoeL, sym(S_n), subs_expr);
                    end
                end
            end
        end

    case 2
        % Unique delta parameters for each interaction type
        for j = 2:level_cut
            for i = 1:count
                V_n = string(base_symvar(i)) + "_" + num2str(j);
                if any(contains(strvar_list, V_n))
                    delta = sym(['delta_', num2str(i)], 'real');
                    coeff = (Rnn(j) - RdL(base_num(i)));
                    subs_expr = base_symvar(i) * exp(-coeff / delta);
                    H_hr(1).HcoeL = subs(H_hr(1).HcoeL, sym(V_n), subs_expr);
                end
            end
        end

        % Substitute overlap terms if present
        if H_hr(1).overlap
            for j = 2:level_cut
                for i = 1:count_S
                    S_n = string(base_symvar_S(i)) + "_" + num2str(j);
                    if any(contains(symvar_list_S, S_n))
                        delta = sym(['delta__2_', num2str(i)], 'real');
                        coeff = (Rnn(j) - RdL(base_num_S(i)));
                        subs_expr = base_symvar_S(i) * exp(-coeff / delta);
                        H_hr(2).HcoeL = subs(H_hr(2).HcoeL, sym(S_n), subs_expr);
                    end
                end
            end
        end
end

% Finalize HR object state
H_hr.num = true;
H_hr.coe = false;
end