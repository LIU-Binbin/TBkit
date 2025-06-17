function SymOper = rotation(angle, axis, inversion, U, spin,options,propArgs)
%ROTATION Generate symmetry operation operator for rotation transformations
%   This function creates a symmetry operator combining spatial rotation and 
%   spin rotation components. Supports both numeric and symbolic calculations.
%
%   Input Arguments:
%       angle       - Rotation angle (in fractions of 2π)
%       axis        - [Optional] 3D rotation axis vector (default: NaN for 2D rotation)
%       inversion   - [Optional] Logical flag for spatial inversion (default: false)
%       U           - [Optional] Predefined rotation matrix or Pauli matrices object
%       spin        - [Optional] Spin quantum number or spin rotation matrix
%       options     - [Optional] Structure with additional parameters:
%           .sym         - Enable symbolic calculation (default: false)
%           .rightorleft - Rotation direction ('left' or 'right', default: 'right')
%       propArgs    - [Optional] Additional properties for Oper class
%
%   addition note for rightorleft
%       right: (1,0) -> (cos(theta),sin(theta))
%   Output:
%       SymOper - Oper object containing combined rotation operation
%
%   Example:
%       % Create a right-handed π/2 rotation about z-axis
%       op = rotation(1/4, [0 0 1]);

arguments
    angle                       % Rotation angle (fraction of 2π)
    axis          = nan         % Rotation axis (NaN for 2D)
    inversion logical = false   % Spatial inversion flag
    U             = nan         % Rotation matrix/Pauli matrices
    spin          = nan         % Spin quantum number/matrix
    options.sym logical = false % Symbolic calculation flag
    options.rightorleft {mustBeMember(options.rightorleft, {'left','right'})} = 'right'
    propArgs.?Oper               % Additional Oper properties
end
propertyCell = namedargs2cell(propArgs);
if isa(U,'double')
    if ~all(all(isnan(U))) && ~all(all(all(isnan(spin))))
        raise ValueError('Only one of `U` and `spin` may be provided.');
    end
elseif isa(U,'pauli_matrix')
    U = double(U);
end
if options.sym
    if isa(angle,'sym')
    else
        angle = sym(2 * pi * angle);
    end
    U_PG = sym(U);
else
    if isa(angle,'sym')
    else
        angle = 2 * pi * angle;
    end
    U_PG = U;
end
if strcmp(options.rightorleft,'left')
    rightorleft = -1;
else
    rightorleft = 1;
end
if isnan(axis)
    R_PG = ([[cos(angle), sin(angle)];
        [-sin(angle), cos(angle)]]);
    if ~isnan(spin)
        U_PG = Oper.spin_rotation(angle * ([0, 0, 1]), spin);
    end
elseif length(axis) == 3
    n = angle * axis / norm(axis);
    R_PG = Oper.spin_rotation(n, rightorleft*Oper.L_matrices(3, 1));
    if inversion
        R_PG = -R_PG;
    else
    end
    if ~isnan(spin)
        if ~isvector(spin)
            U_PG = Oper.spin_rotation(n, spin);
        else
        end
    end
else
    raise ValueError('`axis` needs to be `None` or a 3D vector.')
end
SymOper = Oper(R_PG,U_PG,propertyCell{:});
end
