function str_mat = cell2str(cell_mat)
% CELL2STR Converts a cell array to a string array.
%
% Input:
%   cell_mat - A cell array (m x n) to be converted to a string array.
%
% Output:
%   str_mat - A string array (m x n) corresponding to the input cell array.

    % Apply cellfun to convert each element in the cell array to a string
    str_mat = cellfun(@string, cell_mat, 'UniformOutput', false);
    
    % Convert the result into a string array (m x n)
    str_mat = string(str_mat);
end

% 
% function str_mat = cell2str(cell_mat)
% [m, n] = size(cell_mat);
% str_mat = strings(m,n);
% for i = 1:m
%     for j = 1:n
%         str_mat(i,j) = string(cell_mat{i,j});
%     end
% end
% end