function DOSCAR = DOSCAR_gen(GREENCAR,mode)
%DOSCAR_GEN Generate density of states from Green's functions
%
%   Syntax:
%       DOSCAR = DOSCAR_gen(GREENCAR,mode)
%
%   Description:
%       Computes density of states (DOS) from Green's function data.
%       Supports cell array format for Green's function matrices.
%
%   Inputs:
%       GREENCAR - Green's function data (cell array of matrices)
%       mode     - Calculation mode ('green' for standard DOS)
%
%   Output:
%       DOSCAR - Density of states array
%
%   See also: GREENCAR_gen, EIGENSOLVE
% nargin
if nargin <2
    mode = 'green';
end
if strcmp(mode,'green')
    if iscell(GREENCAR)
        [Nsize1,Nsize2] = size(GREENCAR);
        DOSCAR =zeros(Nsize1,Nsize2);
        for i = 1:Nsize1
            for j = 1:Nsize2
                DOSCAR(i,j) = -trace(imag(GREENCAR{i,j}));
            end
        end
    else
        disp('temp support cell format for GREENCAR');
    end
end
end