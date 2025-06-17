function Am = rotate(A, Rc, Rf, tf, rightorleft, optionsOper, options)
%ROTATE  Rotate a BasisFunc object using specified rotation matrices and translation.
%
%   Am = ROTATE(A, Rc, Rf, tf, rightorleft, optionsOper, options) applies a rotation
%   and translation transformation to the input BasisFunc object array A. The transformation
%   is defined by the following parameters:
%
%       Rc          - A 3x3 rotation matrix for coordinate transformation.
%       Rf          - A 3x3 matrix for transforming additional factors.
%       tf          - A translation vector or matrix (3x1 or 1x3).
%       rightorleft - A string ('right' or 'left') specifying the multiplication order
%                     for the transformation matrices.
%       optionsOper - A structure that may contain an Oper object in its field 'Oper'. If provided,
%                     the Oper object is used to override the parameters Rc, Rf, tf, and associated
%                     flags (conjugate, antisymmetry).
%       options     - A structure with additional options:
%                         sym          - (logical) If true, forces symbolic output (default: false).
%                         conjugate    - (logical) If true, applies conjugation (default: false).
%                         antisymmetry - (logical) If true, applies antisymmetry (default: false).
%                         forgetcoe    - (logical) If true, ignores coefficients (default: false).
%                         fast         - (logical) If true, uses fast computation mode (default: true).
%                         hybird       - (logical) Hybrid mode flag (default: false).
%                         spincoupled  - (logical) If true, accounts for spin coupling (default: false).
%                         orbcoupled   - (logical) If true, accounts for orbital coupling (default: false).
%                         raw          - (logical) If true, returns raw output (default: true).
%                         vpalevel     - (double) Rounding precision for numeric output (default: 6).
%                         center       - (1x3 double) Center of rotation (default: [0,0,0]).
%
%   Inputs:
%       A           - An array of BasisFunc objects.
%       Rc          - 3x3 rotation matrix (default: diag([1, 1, 1])).
%       Rf          - 3x3 matrix (default: diag([1, 1, 1])).
%       tf          - 3x1 or 1x3 translation vector (default: [0, 0, 0]).
%       rightorleft - String specifying transformation order, e.g., 'right' (default).
%       optionsOper - Structure for operator override (default: empty Oper in field Oper).
%       options     - Additional options as described above.
%
%   Output:
%       Am          - The transformed array of BasisFunc objects.
%
%   Example:
%       % Rotate basis functions using default parameters:
%       Am = rotate(A, diag([1, 1, 1]), diag([1, 1, 1]), [0, 0, 0], 'right', struct('Oper', []), struct());
%
%   See also: rotaterow
arguments
    A BasisFunc;
    Rc {mustBeSize(Rc,[3 3])}= diag([1 1 1]);%
    Rf {mustBeSize(Rf,[3 3])}= diag([1 1 1]);%
    tf {mustBeSize(tf,[3 1;1 3])}= ([0 0 0]);%
    rightorleft = 'right';
    optionsOper.Oper = [];
    options.sym = false;
    options.conjugate = false;
    options.antisymmetry = false;
    options.forgetcoe = false;
    options.fast = true;
    options.hybird = false;
    options.spincoupled = false;
    options.orbcoupled = false;
    options.raw = true;
    options.vpalevel = 6;
    options.center = [0,0,0];
end
if ~isempty(optionsOper.Oper)
    Rc = optionsOper.Oper.R;
    Rf = optionsOper.Oper.Rf;
    tf = optionsOper.Oper.tf;
    options.conjugate = optionsOper.Oper.conjugate;
    options.antisymmetry = optionsOper.Oper.antisymmetry;
    optionsCell = namedargs2cell(options);
    Am = rotate(A,Rc,Rf,tf,rightorleft,optionsCell{:});
    return
end
optionsCell = namedargs2cell(options);
Am = rotaterow(A(1,:),Rc,Rf,tf,rightorleft,optionsCell{:});
for i = 2:size(A,1)
    Am = [Am;rotaterow(A(i,:),Rc,Rf,tf,rightorleft,optionsCell{:})];
end
end

