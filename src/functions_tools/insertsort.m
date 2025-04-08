function [sortA, indSortA] = insertsort(A)
    % INSERTSORT - Sorts an array using the insertion sort algorithm.
    % 
    % Parameters:
    %   A : array
    %       The array to be sorted.
    %
    % Returns:
    %   sortA : array
    %       The sorted array.
    %   indSortA : array
    %       The indices that would sort the array A.
    
    len = length(A);  % Length of the input array
    
    % Initialize the index array for sorting
    if isrow(A)
        indSortA = 1:len;  % If A is a row vector, use the same range for indices
    else
        indSortA = (1:len)';  % Ensure the index array is a column vector
    end
    
    % Insertion Sort Algorithm
    for w = 1:len-1
        for v = w+1:-1:2
            try
                % Compare elements in the array
                tmplogical = logical(A(v) < A(v-1)); 
            catch
                % If elements are symbolic, compare them using symbolic sorting
                tmpsym = (A(v) < A(v-1));
                if isa(tmpsym, 'sym')
                    [~, seq_list] = sort([lhs(tmpsym), rhs(tmpsym)]);
                    tmplogical = seq_list(1) < seq_list(2);
                else
                    error('Unknown error in comparison!');
                end
            end
            
            % If the current element is smaller than the previous one, swap them
            if tmplogical
                % Swap elements in A
                tmp = A(v-1);
                tmp_ind = indSortA(v-1);
                A(v-1) = A(v);
                indSortA(v-1) = indSortA(v);
                A(v) = tmp;
                indSortA(v) = tmp_ind;
            else
                % No need to swap, break the loop
                break;
            end
        end
    end
    
    % Return the sorted array and its index order
    sortA = A;
end
