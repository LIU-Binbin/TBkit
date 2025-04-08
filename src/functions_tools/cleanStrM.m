function StrM = cleanStrM(StrM)
% CLEANSTRM Removes all-empty rows and columns from a string matrix.
%
% Input:
%   StrM - A string matrix that may contain empty rows and columns.
%
% Output:
%   StrM - The cleaned string matrix with empty rows and columns removed.

    % Remove empty rows
    cleanRows = ~all(StrM == "", 2); % Check for rows that are not all empty
    StrM = StrM(cleanRows, :);

    % Remove empty columns
    cleanCols = ~all(StrM == "", 1); % Check for columns that are not all empty
    StrM = StrM(:, cleanCols);
end


% function StrM = cleanStrM(StrM)
% 
% cleanlabel = logical(1:size(StrM,1));
% for i = 1:size(StrM,1)
%     if ~isequal(StrM(i,:),repmat("",size(StrM(i,:))))
%         cleanlabel(i) = false;
%     end
% end
% StrM(cleanlabel,:)= [];
% %
% cleanlabel = logical(1:size(StrM,2));
% for j = 1:size(StrM,2)
%     if ~isequal(StrM(:,j),repmat("",size(StrM(:,j))))
%         cleanlabel(j) = false;
%     end
% end
% StrM(:,cleanlabel)= [];
% end