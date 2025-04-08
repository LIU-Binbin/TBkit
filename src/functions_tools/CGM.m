function [Ergebnis, j3, m3] = CGM(j1, j2, m1, m2)
    % Evaluation of all Clebsch-Gordan coefficients for fixed M and
    % m1, m2, j1, j2.
    % Call: [Ergebnis, j3, m3] = CGM(j1, j2, m1, m2)
    % Example: j1 = 1/2, j2 = 1, m3 = 1/2, m1 = 1/2, m2 = 1
    
    % Check if any input values are invalid
    if j1 == -1 || j2 == -1 || isnan(m1) || isnan(m2)
        Ergebnis = nan;
        j3 = -1;
        m3 = nan;
        return;
    end
    
    % Calculate m3 as the sum of m1 and m2
    m3 = m1 + m2;
    
    % Determine the allowed values for j3
    j3 = abs(j1 - j2):j1 + j2;
    j3 = j3(j3 >= abs(m3));   % Only allow j3 values that satisfy |m3| <= j3 <= j1 + j2
    
    % Preallocate the Ergebnis array for efficiency
    Ergebnis = zeros(1, length(j3));
    
    % Loop over the allowed j3 values and compute the Clebsch-Gordan coefficients
    for n = 1:length(j3)
        Ergebnis(1, n) = CGsingle(j1, j2, j3(n), m1, m2, m3);
    end
    
end
