function StrL = char2list(str)
% CHAR2LIST Splits the input string by newlines and removes empty lines.
%
% Input:
%   str - A character array (string) to be split into a list.
%
% Output:
%   StrL - A cell array of strings, with empty lines removed.

    % Split the string by newline characters
    StrL = split(str, '\n');
    
    % Remove empty lines
    StrL = StrL(~cellfun('isempty', StrL));
end


% function StrL = char2list(str)
% arguments
%     str char;
% end
% StrL = split(str,'\n');
% checkL  = logical(1:numel(StrL));
% for i = numel(StrL)
%     if strcmp(StrL(i),'')
%         checkL(i) = false;
%     end
% end
% StrL = StrL(checkL);
% end