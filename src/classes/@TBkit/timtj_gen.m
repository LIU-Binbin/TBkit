function TBkitobj = timtj_gen(TBkitobj, mode)
    % timtj_gen: Generates the tight-binding hopping terms ti - tj for the system.
    %   This function calculates the differences between the orbitals' positions (ti - tj)
    %   in both Cartesian and fractional coordinates, and optionally returns symbolic expressions
    %   for the hopping terms.
    %
    % Inputs:
    %   - TBkitobj: The TBkit object that contains orbital and lattice information.
    %   - mode: A string that determines the output format:
    %           - 'num' (default): Numeric output (Cartesian and fractional).
    %           - 'sym': Symbolic output (for Cartesian and fractional coordinates).
    %
    % Outputs:
    %   - TBkitobj: The updated TBkit object with the calculated hopping terms stored.
    
    if nargin < 2
        mode = 'num';  % Default mode is 'num' for numeric output
    end

    % Initialize the hopping terms (ti - tj) for both Cartesian and fractional coordinates
    % The size of timtj{1} corresponds to the number of orbitals (orbL) and the dimensionality (Dim)
    TBkitobj.timtj{1} = zeros([size(TBkitobj.orbL, 1), size(TBkitobj.orbL, 1), TBkitobj.Dim]);
    
    % If symbolic mode is requested, initialize symbolic arrays
    if strcmp(mode, 'sym')
        TBkitobj.timtj{1} = sym(TBkitobj.timtj{1});
        TBkitobj.timtj{3} = sym(zeros(size(TBkitobj.orbL, 1)));
    end
    
    % Initialize the fractional coordinate version of timtj
    TBkitobj.timtj{2} = TBkitobj.timtj{1};  % Same initial size as timtj{1}
    
    % Calculate the hopping terms (ti - tj) for each pair of orbitals
    for i = 1:size(TBkitobj.orbL, 1)
        ti = TBkitobj.orbL(i, :);  % Get the position of the i-th orbital
        for j = 1:size(TBkitobj.orbL, 1)
            tj = TBkitobj.orbL(j, :);  % Get the position of the j-th orbital
            tij = ti - tj;  % Difference between the i-th and j-th orbital positions
            tij_cart = tij * TBkitobj.Rm;  % Convert the difference to Cartesian coordinates

            % Store the differences in both Cartesian and fractional coordinates
            for d = 1:TBkitobj.Dim
                TBkitobj.timtj{1}(i, j, d) = tij_cart(d);  % Cartesian difference
                TBkitobj.timtj{2}(i, j, d) = tij(d);  % Fractional difference
            end
        end
    end

    % If symbolic mode is requested, calculate the symbolic exponential terms for the hopping terms
    if strcmp(mode, 'sym')
        % Calculate the exponential of the inner product of k and tij in Cartesian coordinates
        ExpInnerTerm = TBkit.matrixtimespage(TBkitobj.VarsSeqLcart(1:TBkitobj.Dim), TBkitobj.timtj{1});
        TBkitobj.timtj{3} = exp(1i * (sum(ExpInnerTerm, 3)));  % Symbolic hopping term (Cartesian)

        % Calculate the exponential of the inner product of k and tij in fractional coordinates
        ExpInnerTermFrac = TBkit.matrixtimespage(TBkitobj.VarsSeqLcart(1:TBkitobj.Dim), TBkitobj.timtj{2});
        TBkitobj.timtj{4} = exp(1i * (sum(ExpInnerTermFrac, 3)));  % Symbolic hopping term (Fractional)
    else
        % For numeric mode, clear symbolic fields
        TBkitobj.timtj{3} = [];
        TBkitobj.timtj{4} = [];
    end
end
