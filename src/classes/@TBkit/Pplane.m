function Pmat = Pplane(A,mode)
            % project a 3D vector to a plane
            if nargin < 2
                mode = 'left';
            end
            if strcmp(mode,'left')
                Pmat = A.'*inv(A*A.')*A;
            else
                Pmat = A*inv(A.'*A)*A.';
            end
        end