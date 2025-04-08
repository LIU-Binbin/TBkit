function [StrL1, StrL2] = ExtractContainPat(StrL, pat, options)
    % EXTRACTCONTAINPAT - Extract strings that contain or do not contain a specified pattern
    %
    % This function takes a vector of strings `StrL` and separates them into two outputs:
    % - `StrL1` contains strings that *do not* contain the specified pattern `pat`.
    % - `StrL2` contains strings that *do* contain the specified pattern `pat`.
    %
    % Input Arguments:
    %   StrL      - A vector of strings to be processed.
    %   pat       - A string pattern to search for within each string (default is "=").
    %   options   - A structure containing optional parameters:
    %       options.IgnoreCase - A logical flag to ignore case during the pattern matching (default is false).
    %
    % Output Arguments:
    %   StrL1     - Strings that do *not* contain the pattern `pat`.
    %   StrL2     - Strings that *contain* the pattern `pat`.
    %
    % Example:
    %   [StrL1, StrL2] = ExtractContainPat(["Hello", "World", "test=123"], "=");
    %   StrL1 will be ["Hello", "World"], and StrL2 will be ["test=123"].
    %
    % Created by: [Your Name]
    
    arguments
        StrL string {mustBeVector};  % Ensure StrL is a vector of strings
        pat = "=";                  % Default pattern to search for is "="
        options.IgnoreCase = false; % Default is case-sensitive matching
    end

    % Initialize a logical array to keep track of which strings contain the pattern
    SelectL = logical(1:length(StrL)); 
    
    % Loop through each string in StrL
    for i = 1:numel(StrL)
        % Check if the current string contains the pattern `pat`
        if contains(StrL(i), pat, "IgnoreCase", options.IgnoreCase)
            % If the pattern is found, mark this index for exclusion in StrL1
            SelectL(i) = false;
        end
    end
    
    % Strings that do *not* contain the pattern (i.e., marked true in SelectL) go into StrL1
    StrL1 = StrL(SelectL);
    
    % Strings that contain the pattern (i.e., marked false in SelectL) go into StrL2
    StrL2 = StrL(~SelectL);
end
