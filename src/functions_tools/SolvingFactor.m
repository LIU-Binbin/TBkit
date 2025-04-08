function plist = SolvingFactor(n)
    % SolvingFactor - Find all factors of a given integer n.
    % 
    % Parameters:
    %   n : integer
    %       The number for which the factors are to be determined.
    %
    % Returns:
    %   plist : array
    %       A list of factors of n, including 1 and n itself.
    
    % Ensure input is a positive integer
    arguments
        n {mustBeInteger, mustBePositive};
    end
    
    % Initialize list of factors
    plist = [1, n];  % 1 and n are always factors
    count = 2;  % Current count of factors found
    p = 2;  % Start checking from 2

    % Loop over potential factors up to sqrt(n)
    while p * p <= n
        if mod(n, p) == 0  % If p is a factor of n
            plist(count + 1) = p;  % Add p to the list
            count = count + 1;
            plist(count + 1) = n / p;  % Add corresponding factor n/p
            count = count + 1;
        end
        p = p + 1;  % Increment p to check the next potential factor
    end
    
    % Remove duplicate factors if any
    plist = unique(plist);
end
