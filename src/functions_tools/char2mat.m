function StrM = char2mat(str)
% CHAR2MAT Converts a character array into a matrix by splitting by newline
% and then processing it with the function `park.str2mat`.
%
% Input:
%   str - A character array to be processed.
%
% Output:
%   StrM - The resulting matrix after processing the split lines.

    % Split the input string by newline character
    StrL = split(str, '\n');
    
    % Call external function to convert the cell array to a matrix
    StrM = str2mat_overwrited(StrL);
end
