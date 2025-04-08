function [parity] = Parity(filename, Occupiband)
    % Parity computes the total, positive, and negative contributions
    % from a data file, adjusting for a specified occupancy band.
    %
    % Inputs:
    %   filename   - Name of the text file containing the data.
    %   Occupiband - Optional parameter specifying the occupancy band.
    %                If not provided, the function attempts to read 'Efermi' from a file.
    %
    % Outputs:
    %   parity     - Structure containing the total, positive, and negative contributions.
    
    % Set default value for Occupiband if not provided
    if nargin < 2
        % Attempt to read 'Efermi' from a file
        Efermi = textread('Efermi');
        % Compute geigval and find the occupancy sequence
        geigval = double(dataArray{3}) - Efermi;
        Occupiseq = find(geigval > 0, 1) - 1;
    else
        % Read data from the specified file
        data = readtable(filename, 'Format', '%3s%3s%11s%11s%12s');
        % Convert relevant columns to numeric values
        bnd = double(data.Var1);
        I = double(data.Var5);
        % Find the occupancy sequence based on the provided Occupiband
        Occupiseq = find(bnd > Occupiband, 1) - 1;
    end
    
    % Initialize variables
    parity.tot = sum(I(1:Occupiseq));
    Minus = 0;
    Plus = 0;
    
    % Compute positive and negative contributions
    for i = 1:Occupiseq
        if I(i) > 0
            Plus = Plus + I(i);
        elseif I(i) < 0
            Minus = Minus + I(i);
        end
    end
    
    % Adjust contributions if necessary
    if -Minus + Plus < Occupiband
        temp_num = (Occupiband + Minus - Plus) / 2;
        Minus = Minus - temp_num;
        Plus = Plus + temp_num;
    end
    
    % Assign results to output structure
    parity.minus = Minus;
    parity.plus = Plus;
end
