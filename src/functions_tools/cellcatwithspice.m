function OutputStr = cellcatwithspice(StrCell)
% CELLCATWITHSPICE Concatenate strings from a cell array into a single string
% with spaces between elements and newline at the end of each row.
%
% Input:
%   StrCell - A cell array containing strings or a single string.
%
% Output:
%   OutputStr - A concatenated string with elements from StrCell, separated
%               by spaces, and newlines after each row.

    % If input is a single string, return it with a newline
    if ischar(StrCell)
        OutputStr = strcat(StrCell, "\n");
        return;
    end

    % Initialize the output string
    OutputStr = '';

    % Loop through each row in the cell array
    for i = 1:length(StrCell)
        % Concatenate the current row's elements with spaces between them
        rowStr = strjoin(StrCell{i}, ' ');
        
        % Append the row string followed by a newline to the output
        OutputStr = strcat(OutputStr, rowStr, "\n");
    end
end

% function OutputStr = cellcatwithspice(StrCell)
% OutputStr = '';
% if ischar(StrCell)
%     OutputStr = strcat(OutputStr,StrCell,"\n");
%     return;
% end
% for i = 1:length(StrCell)
%     for j = 1:length(StrCell{i})-1
%         OutputStr = strcat(OutputStr,StrCell{i}(j)," ");
%     end
%     OutputStr = strcat(OutputStr,StrCell{i}(j+1),"\n");
% end
% end