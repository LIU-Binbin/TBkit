function Pmat = Pvector(A)
            if nargin < 2
                mode = 'left';
            end
            if strcmp(mode,'left')
                Pmat = A.'*A/(A*A.');
            else
                Pmat = A*A.'/(A.'*A);
            end
        end