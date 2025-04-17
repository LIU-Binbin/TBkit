function [list_obj_unique, sorted_label, cut_slice] = cut_tools(list_obj)
%CUT_TOOLS Sort, deduplicate, and identify index ranges for repeated rows.
%
%   [list_obj_unique, sorted_label, cut_slice] = cut_tools(list_obj)
%
%   Inputs:
%       list_obj       - NxM array, where each row is a vector (e.g., spatial coordinate)
%
%   Outputs:
%       list_obj_unique - KxM array, unique rows from sorted list
%       sorted_label    - Nx1 vector, index mapping from original to sorted list
%       cut_slice       - Kx2 array, start and end index in sorted list for each unique row

    % Sort the rows of input matrix
    [list_obj_sorted, sorted_label] = sortrows(list_obj);

    % Find unique rows and their first occurrence in the sorted list
    [list_obj_unique, unique_label] = unique(list_obj_sorted, 'rows', 'stable');

    % Determine index range for each group of identical rows in the sorted list
    group_start = unique_label;
    group_end = [unique_label(2:end)-1; size(list_obj_sorted,1)];
    cut_slice = [group_start, group_end];
end
