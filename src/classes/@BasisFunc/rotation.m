function [U,orbL] = rotation(A, Rc, Rf, tf, optionsConvection, optionsOper, optionsRm, options,optionsAddition)
%ROTATION  Rotate a BasisFunc object and compute the inner product with the original.
%
%   U = ROTATION(A, Rc, Rf, tf, optionsConvection, optionsOper, optionsRm, options)
%   applies a series of transformations to the input BasisFunc object A by using 
%   the specified rotation matrices (Rc for coordinate rotation, Rf for orbital rotation), 
%   a translation vector tf, and additional options controlling convection, operator 
%   overrides, and removal functionality. After rotation, the function computes the inner 
%   product of the rotated object with the original A, yielding the transformation matrix U.
%
%   Inputs:
%       A                   - A BasisFunc object.
%       Rc                  - A 3x3 rotation matrix for coordinate transformation 
%                             (default: diag([1 1 1])).
%       Rf                  - A 3x3 rotation matrix for orbital transformation 
%                             (default: diag([1 1 1])).
%       tf                  - A translation vector (or matrix) (default: [0 0 0]).
%       optionsConvection   - A structure containing:
%                               rightorleft: 'right' or 'left' (default: 'right').
%       optionsOper         - A structure for operator overrides with field 'Oper'
%                             (default: empty Oper object).
%       optionsRm           - A structure with field Rm used for removal operations.
%                             (default: POSCAR_read).
%       options             - A structure with additional options:
%                               sym           - (logical) Force symbolic computations (default: false).
%                               conjugate     - (logical) Apply conjugation (default: false).
%                               antisymmetry  - (logical) Apply antisymmetry (default: false).
%                               forgetcoe     - (logical) Ignore coefficients (default: false).
%                               fast          - (logical) Use fast computation mode (default: true).
%                               hybird        - (logical) Use hybrid mode (default: false).
%                               spincoupled   - (logical) Consider spin coupling (default: false).
%                               orbcoupled    - (logical) Consider orbital coupling (default: false).
%                               raw           - (logical) Return raw output (default: true).
%                               vpalevel      - (double) Rounding precision (default: 6).
%                               center        - (1x3 double) Center of rotation (default: [0,0,0]).
%
%   Output:
%       U  - Transformation matrix computed as the inner product between the rotated 
%            BasisFunc object and the original A.
%
%   Example:
%       % Rotate the BasisFunc object A using specified rotation matrices and compute U:
%       optsConv.rightorleft = 'right';
%       optsOper.Oper = [];        % No operator override
%       optsRm.Rm = POSCAR_read;    % Default removal matrix from POSCAR_read
%       opts.sym = false; opts.conjugate = false; opts.antisymmetry = false;
%       U = rotation(A, diag([1 1 1]), diag([1 1 1]), [0 0 0], optsConv, optsOper, optsRm, opts);
%
%   See also: rotate, InnerProduct, BasisFunc

    arguments
        A BasisFunc;
        Rc {mustBeSize(Rc, [3 3])} = diag([1 1 1]);
        Rf {mustBeSize(Rf, [3 3])} = diag([1 1 1]);
        tf {mustBeSize(tf, [3 1;1 3])} = ([0 0 0]);
        optionsConvection.rightorleft {mustBeMember(optionsConvection.rightorleft, {'right','left'})} = 'right';
        optionsOper.Oper = [];
        optionsRm.Rm = POSCAR_read;
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
        optionsAddition.IgnoreOrbL = false;
    end

    rightorleft = optionsConvection.rightorleft;
    optionsCell = namedargs2cell(options);
    if ~isempty(optionsOper.Oper)
        Rf = optionsOper.Oper.Rf;
        if isempty(Rf)
            optionsOper.Oper = optionsOper.Oper.attachRm(optionsRm.Rm);
        end
        Rf = optionsOper.Oper.Rf;
        Rc = optionsOper.Oper.R;
        tf = optionsOper.Oper.tf;% seize notation; t = tf
        options.conjugate = optionsOper.Oper.conjugate;
        options.antisymmetry = optionsOper.Oper.antisymmetry;
        optionsCell = namedargs2cell(options);
        Am = rotate(A, Rc, Rf, tf, rightorleft, optionsCell{:}); % Oper->Basis
    else
        Am = rotate(A, Rc, Rf, tf, rightorleft, 'Oper', optionsOper.Oper, optionsCell{:});
    end
    if optionsAddition.IgnoreOrbL
       for i = 1:length(A)
           A(i).BForb = Am(i).BForb;
           orbL(i,:) = Am(i).BForb;
       end
    else
        orbL = [];
    end
    U = InnerProduct(Am, A, 'sym', options.sym);% <A,Oper,A> 
end

