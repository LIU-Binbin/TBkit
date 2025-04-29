function [c,t] = coeffsAll(p,varargin)
            [Ni,Nj] = size(p);
            for i = 1:Ni
                for j = 1:Nj
                    switch length(varargin)
                        case 0
                            [c{i,j},t{i,j}] = coeffs(p(i,j));
                        otherwise
                            [c{i,j},t{i,j}] = coeffs(p(i,j),varargin);
                    end
                end
            end
        end