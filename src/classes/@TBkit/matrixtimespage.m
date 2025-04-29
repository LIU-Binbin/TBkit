function pagenew = matrixtimespage(mat,page)
            if length(mat) ~= size(page,3)
                raise ValueError('the length of the mat must equals with the pages.')
            end
            pagenew = page;
            if isvector(mat)
                for i = 1:length(mat)
                    pagenew(:,:,i) = mat(i)*page(:,:,i);
                end
            else
                for i = 1:length(mat)
                    pagenew(:,:,i) = sum(TBkit.matrixtimespage(mat(i,:),page),3);
                end
            end
        end