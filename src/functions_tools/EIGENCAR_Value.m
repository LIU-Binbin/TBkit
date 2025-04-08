function ValueTotal = EIGENCAR_Value(EIGENCAR_DFT, EIGENCAR, extra_parm, options)
    % EIGENCAR_Value - Computes the comparison value between two eigenvalue matrices: EIGENCAR_DFT and EIGENCAR.
    %
    % This function calculates the difference between two eigenvalue datasets using various metrics
    % (such as direct differences and differences of differences). The results are weighted by the provided 
    % 'extra_parm' vector and summed up to produce a final total comparison value.
    %
    % Input Arguments:
    %   EIGENCAR_DFT      - Eigenvalue matrix from the DFT calculation or reference data (double).
    %   EIGENCAR          - Eigenvalue matrix from another method or model (double).
    %   extra_parm        - Vector of scaling factors for the computed metrics. Default is [1,1] (double).
    %   options           - Structure containing options for the comparison. Fields:
    %       options.mode     - Comparison mode, can be 'default' or 'extra'. Default is 'default'.
    %       options.algorithm - The algorithm used for comparison. Default is 'pure_comparision'.
    %
    % Output:
    %   ValueTotal        - The total comparison value computed from weighted differences.
    %
    % Example:
    %   ValueTotal = EIGENCAR_Value(EIGENCAR_DFT, EIGENCAR, [1,2], options);
    
    arguments
        EIGENCAR_DFT double;   % Eigenvalues from the DFT (reference) calculation
        EIGENCAR double;       % Eigenvalues from the model (comparison)
        extra_parm double = [1,1];  % Weighting parameters for the metrics (default is [1,1])
        options.mode = 'default';  % Mode for comparison: 'default' or 'extra'
        options.algorithm = 'pure_comparision';  % Algorithm for comparison (default is 'pure_comparision')
    end

    % If 'extra' mode is selected, retrieve additional parameters from the base workspace
    if strcmp(options.mode, 'extra')
        options_extra = evalin('base', 'options_extra');  % Fetch additional parameters from the base workspace
    end
    
    % Select the data for comparison based on the chosen algorithm and mode
    if strcmp(options.algorithm, 'pure_comparision') && ~strcmp(options.mode, 'extra')
        DATA1 = EIGENCAR_DFT;  % Reference eigenvalue data (DFT results)
        DATA2 = EIGENCAR;      % Eigenvalue data from the model (to compare)
    elseif strcmp(options.algorithm, 'pure_comparision') && strcmp(options.mode, 'extra')
        % Use extra parameters (e.g., range of bands and k-points) for comparison
        NBAND_range = options_extra.NBAND_range;  % Band range from extra options
        klist_range = options_extra.klist_range;  % k-point range from extra options
        DATA1 = EIGENCAR_DFT(NBAND_range, klist_range);  % Subset of the reference data
        DATA2 = EIGENCAR(NBAND_range, klist_range);      % Subset of the comparison data
    end
    
    % Initialize an array to store the comparison values (one for each metric)
    metrics = zeros(1, 2);  % Pre-allocate array for two metrics
    
    % Direct difference comparison: root mean squared error (RMSE)
    % This computes the square root of the mean squared difference between the two datasets
    metrics(1) = extra_parm(1) * sqrt(mean(mean(abs(DATA1 - DATA2)).^2));
    
    % Difference of differences: RMSE of the differences between consecutive eigenvalues
    % This computes the difference of differences and then calculates the RMSE
    metrics(2) = extra_parm(2) * sqrt(mean(mean(abs(diff(DATA1.') - diff(DATA2.')))));
    
    % The final total comparison value is the sum of the weighted metrics
    ValueTotal = sum(metrics);
end
