function A_Lj = rotaterow(A, Rc, Rf, tf, rightorleft, options)
%ROTATEROW  Rotate a row of a BasisFunc object with specified rotation matrices and translation.
%
%   A_Lj = ROTATEROW(A, Rc, Rf, tf, rightorleft, options) applies rotation and translation
%   transformations to a BasisFunc object (or row thereof) using the provided rotation matrices
%   for coordinates (Rc) and for orbitals (Rf), a translation vector (tf), and a specified order
%   ('right' or 'left'). The function supports additional options controlling symbolic behavior,
%   conjugation, antisymmetry, and whether to process spin and/or orbital coupling.
%
%   Inputs:
%       A           - A BasisFunc object (or array) representing basis functions.
%       Rc          - A 3x3 rotation matrix for the coordinate variables (default: diag([1 1 1])).
%       Rf          - A 3x3 rotation matrix for orbital variables (default: diag([1 1 1])).
%       tf          - A translation vector or matrix (3x1 or 1x3) (default: [0 0 0]).
%       rightorleft - A string specifying the multiplication order ('right' is default).
%       options     - A structure with optional fields:
%                        sym           - (logical) Force symbolic computations (default: true).
%                        conjugate     - (logical) Apply conjugation (default: false).
%                        antisymmetry  - (logical) Apply antisymmetry (default: false).
%                        forgetcoe     - (logical) Ignore coefficients (default: false).
%                        fast          - (logical) Use fast computation mode (default: true).
%                        hybird        - (logical) Use hybrid mode (default: false).
%                        spincoupled   - (logical) Consider spin coupling (default: false).
%                        orbcoupled    - (logical) Consider orbital coupling (default: false).
%                        raw           - (logical) Return raw output (default: true).
%                        vpalevel      - (double) Rounding precision for numerical output (default: 6).
%                        center        - (1x3 double) Center of rotation (default: [0,0,0]).
%
%   Output:
%       A_Lj        - The rotated BasisFunc object.
%
%   Example:
%       % Rotate a BasisFunc object A with default parameters:
%       A_Lj = rotaterow(A, diag([1,1,1]), diag([1,1,1]), [0,0,0], 'right', struct());
%
%   See also: rotate, contractrow, BasisFunc.rotation_orb, rotaterow (recursive call)

arguments
    A BasisFunc;
    Rc {Spin.mustBeSize(Rc, [3 3])} = diag([1 1 1]);
    Rf {Spin.mustBeSize(Rf, [3 3])} = diag([1 1 1]);
    tf {Spin.mustBeSize(tf, [3 1;1 3])} = ([0 0 0]);
    rightorleft = 'right';
    options.sym = true;
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

optionsCell = namedargs2cell(options);
if options.hybird
    error('not support yet');
    return
else
    BFuncLtmp = ([A.BFuncL]);
    coeLtmp = ([A.coe]);
    BFuncLtmp = BasisFunc.introduce_coe(BFuncLtmp, coeLtmp);
end

if iscell(A(1).BFuncL)
    error('not support yet');
    return;
end

if ~options.orbcoupled
    orbL = BasisFunc.rotation_orb(A(1).BForb, Rf.', tf, optionsCell{:});
    if ~options.spincoupled
        if isa(A(1).BFuncL, 'Qnum')
            BFuncLtmp = rotaterow(BFuncLtmp, Rc, rightorleft, ...
                'sym', options.sym, 'antisymmetry', options.antisymmetry, 'conjugate', options.conjugate);
            BFuncLtmp = contractrow(BFuncLtmp);
            spinL = Spin([BFuncLtmp.s], [BFuncLtmp.sz]);
            [coeLtmp, BFuncLtmp] = BasisFunc.extract_coe(BFuncLtmp, 'sym', options.sym, 'vpalevel', options.vpalevel);
            A_Lj = BasisFunc(BFuncLtmp, spinL, 1, coeLtmp, orbL);
        else
            BFuncLtmp = rotaterow(BFuncLtmp, Rc, rightorleft, ...
                'sym', options.sym, 'antisymmetry', options.antis
