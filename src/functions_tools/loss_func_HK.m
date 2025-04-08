function Value = loss_func_HK(parameter, extra_param, options)
% loss_func_HK calculates the loss value based on different parameter types and methods.
% It supports three modes of operation: 'Single', 'NM', and 'Bayes'.
% The function uses the Hamiltonian matrix and DFT eigenvalues to compute the loss.

arguments
    parameter;                          % Input parameters for loss calculation
    extra_param double = 0;             % Extra parameter for the loss function (default: 0)
    options.mode char = 'default';      % Mode of operation (default: 'default')
    options.algorithm char = 'pure_comparision'; % Algorithm to use for comparison (default: 'pure_comparision')
end

% Retrieve data from the base workspace
EIGENCAR_DFT = evalin('base', 'EIGENCAR_DFT');
H_hk_n = evalin('base', 'H_hk_n');

% Determine the mode of the parameter
if isa(parameter, 'struct')
    mode = 'Bayes';  % Bayesian mode
elseif isa(parameter, 'double')
    mode = 'NM';     % Numerical minimization mode
    if length(parameter) <= 1
        mode = 'Single';  % Single parameter mode
    end
end

% Main switch for different modes
switch mode
    case 'Single'
        % No specific code for 'Single' mode yet, can add logic if needed.
        
    case 'NM'
        % In Numerical Minimization (NM) mode, update H_hk_n with parameter values
        Varlist = evalin('base', 'Varlist');
        % Substitute the parameters into the Hamiltonian
        H_hk_n = H_hk_n.subs(Varlist, parameter);
        % Perform substitution of all variables
        H_hk_n = H_hk_n.Subsall();
        % Generate eigenvalues for the Hamiltonian
        EIGENCAR_HK = H_hk_n.EIGENCAR_gen();
        % Calculate loss by comparing the generated eigenvalues with the DFT values
        Value = EIGENCAR_Value(EIGENCAR_DFT, EIGENCAR_HK, extra_param, 'mode', options.mode, 'algorithm', options.algorithm);
        
    case 'Bayes'
        % In Bayesian mode, substitute the parameters dynamically
        Varlist = evalin('base', 'Varlist');
        Varlock = evalin('base', 'Varlock');
        
        % Loop through the Varlist to update H_hk_n
        for i = 1:length(Varlist)
            try
                % Try to substitute the parameter from the input structure
                H_hk_n = H_hk_n.subs(Varlist(i), parameter.(string(Varlist(i))));
            catch
                % If the field does not exist in the input, use the locked value
                H_hk_n = H_hk_n.subs(Varlist(i), Varlock(i));
            end
        end
        % Perform substitution of all variables
        H_hk_n = H_hk_n.Subsall();
        % Generate eigenvalues for the Hamiltonian
        EIGENCAR_HK = H_hk_n.EIGENCAR_gen();
        % Calculate loss by comparing the generated eigenvalues with the DFT values
        Value = EIGENCAR_Value(EIGENCAR_DFT, EIGENCAR_HK, extra_param, 'mode', options.mode, 'algorithm', options.algorithm);
end

end
