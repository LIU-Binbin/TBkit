function Parity_list = Parity_eight(Occupi, list)
    % Generate a list of parity values for different label files.
    % Inputs:
    %   Occupi - Occupation states or any input required by the Parity function.
    %   list   - A list of labels (default is [1, 2, ..., 8]).
    %
    % Output:
    %   Parity_list - A structure with parity calculations for each label.
    
    % Set default value for 'list' if not provided
    if nargin < 2
        list = 1:8;  % Default list of labels from 1 to 8
    end
    
    % Initialize an empty structure to store the parity results
    Parity_list = struct();
    
    % Loop through the list of labels
    for i = 1:length(list)
        % Generate the label string (e.g., 'label_file_1', 'label_file_2', etc.)
        label = sprintf('label_file_%d', list(i));
        
        % Store the result in the structure
        Parity_list.(label) = Parity(label, Occupi);
    end
end
