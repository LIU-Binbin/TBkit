function sym_out = double2sym(double_in, accuracy)
    % DOUBLE2SYM Convert double-precision number to symbolic expression.
    %   SYM_OUT = DOUBLE2SYM(DOUBLE_IN) converts a double-precision number to a
    %   symbolic expression using default accuracy (5 decimal places).
    %
    %   SYM_OUT = DOUBLE2SYM(DOUBLE_IN, ACCURACY) specifies the number of decimal
    %   places to consider for rounding.
    %
    %   Example:
    %       double_in = 0.333333333333;
    %       sym_out = double2sym(double_in, 3);
    %
    %   See also: DOUBLE2SYM_ONE

    if nargin < 2
        accuracy = 5;  % Default accuracy if not specified
    end

    if isnan(double_in)
        sym_out = nan; % Return NaN if input is NaN
        return
    end

    % Calculate different powers of the input number
    double_single = roundn(double_in^1, -accuracy);
    double_square = roundn(double_in^2, -accuracy);
    double_cubic = roundn(double_in^3, -accuracy);

    % Convert each power to symbolic form
    double_single_out = double2sym_one(double_single, accuracy);
    double_square_out = double2sym_one(double_square, accuracy);
    double_cubic_out = double2sym_one(double_cubic, accuracy);

    % Create an array of symbolic expressions
    double_out = [double_single_out; double_square_out; double_cubic_out];

    % Calculate the lengths of each symbolic expression
    double_out_length = [...
        length(char(string(double_single_out))); ...
        length(char(string(double_square_out))); ...
        length(char(string(double_cubic_out)))];

    % Find the shortest symbolic expression
    [~, min_label] = min(double_out_length);

    % Determine the final symbolic output
    sym_out = double_out(min_label)^(1/min_label);
end

function sym_out = double2sym_one(double_in, accuracy)
    % DOUBLE2SYM_ONE Helper function to convert a double to a symbolic expression.
    %   SYM_OUT = DOUBLE2SYM_ONE(DOUBLE_IN, ACCURACY) converts a double-precision
    %   number to a symbolic expression, considering the specified accuracy.
    
    integer = fix(double_in);  % Extract integer part
    digit = double_in - integer;  % Extract fractional part

    % Create symbolic representation
    sym_out = sym(integer) + 1 / sym(roundn(1 / digit, -accuracy + 2));
end

