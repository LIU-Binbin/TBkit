function Ergebnis = CGsingle(j1, j2, j3, m1, m2, m3)
    % Evaluation of Clebsch-Gordan coefficient based on all angular values
    % CG = CGsingle(j1,j2,j3,m1,m2,m3)
    
    % Based on Vers. 2.2 18.08.2016
    % Wolfgang Schweizer
    % Simulation of physical systems - Computational Physics with MATLAB
    % Version 1.1 18.10.2019

    Ergebnis = [];

    % Precompute some factors related to j3
    vorfakj3 = [j1+j2-j3, j3+j1-j2, j3+j2-j1];
    vorfakj3 = factorial(vorfakj3);
    vorfakj3 = sqrt((2 * j3 + 1) * prod(vorfakj3) / factorial(j3 + j1 + j2 + 1));

    % Regge symbol check
    regge = [-j1 + j2 + j3, j1 - j2 + j3, j1 + j2 - j3; ...
             j1 - m1, j2 - m2, j3 + m3; ...
             j1 + m1, j2 + m2, j3 - m3];
         
    % Check if any of the values in regge are negative or do not satisfy the sum conditions
    if sum(regge(:) < 0) || any(sum(regge, 2) ~= (j1 + j2 + j3))
        disp('Null!! Incorrect values');
        Ergebnis = [Ergebnis; 0];
        return;
    else
        % Start of the loop to calculate Clebsch-Gordan coefficient
        lauf1 = [j1 + j2 - j3, j1 - m1, j2 + m2, j3 - j2 + m1, j3 - j1 - m2];
        s = abs(min(lauf1));  % Get the minimum of lauf1 and take its absolute value
        laufwhile = true;      % Initialize while loop condition
        notstop = 1;           % Counter for the loop
        vorfak = [j1 + m1, j1 - m1, j2 + m2, j2 - m2, j3 + m3, j3 - m3];
        
        % Check for negative factors in vorfak
        if any(vorfak < 0)
            disp('Erroneous Prefactor');
        else
            % Compute factorial for prefactors
            vorfak = factorial(vorfak);
            vorfak = exp(sum(log(vorfak) * 0.5));  % Numerically stable product calculation
        end
        
        % Initialize an array for storing sum results
        sums = zeros(1, 100);  % Preallocate the sums array
        
        % Begin the while loop for computing Clebsch-Gordan coefficient
        while laufwhile
            test = [j1 + j2 - j3 - s, j1 - m1 - s, j2 + m2 - s, ...
                    j3 - j2 + m1 + s, j3 - j1 - m2 + s];
            testok = min(test);  % Check the minimum value in the test array
            zwi1 = [s, test];
            zwi1 = factorial(zwi1);  % Compute the factorial of each term in zwi1
            zwi1 = exp(sum(log(zwi1)));  % Numerically stable product calculation
            
            % Store the result in the sums array
            sums(notstop) = (-1)^s * vorfak / zwi1;
            s = s + 1;  % Increment s for the next iteration
            
            % Update the test condition
            test = [j1 + j2 - j3 - s, j1 - m1 - s, j2 + m2 - s, ...
                    j3 - j2 + m1 + s, j3 - j1 - m2 + s];
            testok = min(test);
            notstop = notstop + 1;
            
            % Prevent infinite loops by limiting the number of iterations
            if notstop > 100
                disp('While-loop ran 100 times - error?');
                return;
            end
            
            laufwhile = testok >= 0;  % Continue loop if testok is non-negative
        end
        
        % Sum the results of the sums array
        Summe = sum(sums);
        Ergebnis = [Ergebnis; vorfakj3 * sym(Summe)];  % Store the final result
    end
end
