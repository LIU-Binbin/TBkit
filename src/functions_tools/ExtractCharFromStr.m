function [CharL, StrL] = ExtractCharFromStr(StrL, number)
    % EXTRACTCHARFROMSTR - Extracts a specific number of characters from each string in a list
    %
    % This function takes a string array `StrL` and extracts the first `number` characters 
    % from each string in `StrL`. The extracted characters are returned in `CharL`, 
    % and the remaining characters from each string are returned in `StrL`.
    %
    % Input Arguments:
    %   StrL    - A string array (vector of strings).
    %   number  - The number of characters to extract from each string (default is 1).
    %
    % Output:
    %   CharL   - A character array containing the extracted characters.
    %   StrL    - The updated string array with the extracted characters removed.
    %
    % Example:
    %   [CharL, StrL] = ExtractCharFromStr(["Hello", "World"], 2);
    %   CharL will be 'He' and 'Wo', while StrL will be 'llo' and 'rld'.

    arguments
        StrL string {mustBeVector};  % Ensure StrL is a vector of strings
        number = 1;                  % Default value for 'number' is 1 if not provided
    end
    
    % Convert the string array to a cell array of character arrays
    CharCell = convertStringsToChars(StrL);

    % If CharCell is a single character array (non-cell input), extract directly
    if ischar(CharCell)
        % Extract 'number' characters from the start
        CharL = CharCell(1:number);
        % Remove the extracted characters from CharCell
        CharCell(1:number) = [];
        % Convert the remaining characters back to a string
        StrL = convertCharsToStrings(CharCell);
        return;
    end

    % Initialize CharL as a character array, sized to match the input strings
    CharL = repmat("", [numel(CharCell), number]);
    
    % Loop over each string in CharCell and extract the required number of characters
    for i = 1:numel(CharCell)
        CharL(i) = CharCell{i}(1:number);     % Extract characters
        CharCell{i}(1:number) = [];           % Remove extracted characters from the string
    end

    % Convert CharL to a character array
    CharL = char(CharL);

    % Convert the remaining characters in CharCell back to a string array
    StrL = convertCharsToStrings(CharCell);
end
