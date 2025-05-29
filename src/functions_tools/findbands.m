function [allk_min, anyk_max] = findbands(EIGENCAR, Erange)
    % FINDBANDS - Identifies the bands within a specified energy range and 
    % finds the k-point with the minimum and maximum number of bands in that range.
    %
    % Syntax:
    %   [allk_min, anyk_max] = findbands(EIGENCAR, Erange)
    %
    % Inputs:
    %   EIGENCAR - A matrix of eigenvalues (bands x k-points).
    %   Erange   - Energy range [Emin, Emax] to find bands. If not provided,
    %              the user will be prompted to input the values.
    %
    % Outputs:
    %   allk_min - A struct containing the k-point with the minimum number of
    %              bands in the range, and its associated band indices.
    %   anyk_max - A struct containing the k-point with the maximum number of
    %              bands in the range, and its associated band indices.
    arguments
        EIGENCAR 
        Erange 
    end
    % Get the number of bands (Nbands) and k-points (Ktotals)
    [Nbands, Ktotals] = size(EIGENCAR);
    fprintf("%d bands in %d Kpoints ; (Efermi set to zero)\n\n ", Nbands, Ktotals);
    
    % Prompt for energy range if not provided
    if nargin < 2
        Erange(1) = input('Please input Emin (eV): ');
        Erange(2) = input('Please input Emax (eV): ');
    end

    Emin = Erange(1);
    Emax = Erange(2);
    
    % Initialize the findbands matrix to store band indices
    findbands = zeros(Ktotals, 3);
    
    % Loop over each k-point to find the bands in the energy range
    for i = 1:Ktotals
        Minflag = 1;
        Maxflag = 1;
        for j = 1:Nbands
            if EIGENCAR(j, i) > Emin && Minflag == 1
                Minflag = 0;
                findbands(i, 1) = j; % Store the first band index
            end
            if EIGENCAR(j, i) > Emax && Maxflag == 1
                Maxflag = 0;
                findbands(i, 2) = j - 1; % Store the last band index
                findbands(i, 3) = j - findbands(i, 1); % Number of bands in range
                break; % Exit the loop as we've found the range
            end
        end
        % If no Maxband found, set the last band index to the total number of bands
        if Maxflag == 1
            findbands(i, 2) = Nbands;
            findbands(i, 3) = Nbands - findbands(i, 1);
        end
    end

    % Display the findbands matrix if Erange was not provided
    if nargin < 2
        disp(findbands);
    end
    
    % Find the k-point with the minimum number of bands in the range
    [allk_min.nbands, allk_min.label] = min(findbands(:, 3));
    allk_min.n1 = findbands(allk_min.label, 1); % First band index
    allk_min.n2 = findbands(allk_min.label, 2); % Last band index
    
    % Find the k-point with the maximum number of bands in the range
    [anyk_max.nbands, anyk_max.label] = max(findbands(:, 3));
    anyk_max.n1 = findbands(anyk_max.label, 1); % First band index
    anyk_max.n2 = findbands(anyk_max.label, 2); % Last band index
    
    % Print the results
    fprintf("For all k: %d %d %d \n", allk_min.n1, allk_min.n2, allk_min.nbands);
    fprintf("For any k: %d %d %d \n", anyk_max.n1, anyk_max.n2, anyk_max.nbands);
end


