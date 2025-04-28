function pagenew = matrix_mtimes_page(mat,page)
            if length(mat) ~= size(page,2)
                raise ValueError('the length of the mat must equals with the pages.')
            end
            pagenew = page;
            for i = 1:size(page,3)
                pagenew(:,:,i) = mat*page(:,:,i);
            end
        end