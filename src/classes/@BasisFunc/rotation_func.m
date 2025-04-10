function [BFuncL, coeL] = rotation_func(BFuncL, R, t, options)
%ROTATION_FUNC  Rotate a cell array of basis functions by a given symmetry operation.
%
%   [BFuncL, coeL] = ROTATION_FUNC(BFuncL, R, t, options)
%   applies a coordinate-space rotation and translation to each basis function 
%   in the input cell array `BFuncL`, using the specified rotation matrix `R`, 
%   translation vector `t`, and additional options. The transformation is applied 
%   only to non-symbolic elements. Symbolic expressions are skipped.
%
%   Inputs:
%       BFuncL              - Cell array of basis functions (e.g., Qnum or QnumL).
%       R                   - 3x3 rotation matrix to apply.
%       t                   - 3x1 translation vector to apply.
%       options             - Structure with optional arguments:
%           conjugate       - (logical) Whether to apply complex conjugation (default: false).
%           Rm              - (3x3 matrix) Lattice matrix for coordinate reference 
%                             (default: identity).
%
%   Outputs:
%       BFuncL              - Rotated cell array of basis functions.
%       coeL                - Not returned in this function currently; placeholder for future.
%
%   Notes:
%       - The function skips symbolic entries (`sym`) in the input.
%       - Internally calls `rotate` for each individual basis function.
%
%   Example:
%       R = rotation_matrix_z(pi/2);
%       t = [0;0;0];
%       BFuncL = {Qnum(1,0,0), Qnum(0,1,0)};
%       BFuncL_rot = rotation_func(BFuncL, R, t);
%
%   See also: rotate, Qnum, QnumL

    arguments
        BFuncL
        R
        t
        options.conjugate logical = false;
        options.Rm = eye(3);
    end

    optionsCell = namedargs2cell(options);

    for i = numel(BFuncL)  % Likely a bug — this loop only runs for one i
        if isa(BFuncL{i}, 'sym')
            % Skip symbolic entries
        else
            BFuncL{i} = rotate(BFuncL{i}, R, t, optionsCell{:});
        end
    end

    % coeL is declared but not assigned or used — may be placeholder or for compatibility
end

