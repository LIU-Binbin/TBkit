function B = contractrow(A, options)
%CONTRACTROW  Contract rows of a BasisFunc object.
%
%   B = CONTRACTROW(A) processes the rows of the BasisFunc object A by
%   cleaning and contracting its basis function list, while introducing
%   coefficients into the basis. If A contains cell arrays in its BFuncL
%   property or if hybrid mode is enabled, the function errors.
%
%   B = CONTRACTROW(A, options) allows specification of additional processing
%   options. The available options are:
%       forgetcoe     - Logical flag to ignore coefficients (default: false).
%       fast          - Logical flag to use fast computation (default: true).
%       hybird        - Logical flag for hybrid contraction (default: false).
%       spincoupled   - Logical flag for spin coupling (default: false).
%       orbcoupled    - Logical flag for orbital coupling (default: false).
%       raw           - Logical flag to return raw result (default: true).
%       sym           - Logical flag to force symbolic output (default: false).
%       conjugate     - Logical flag to apply conjugation (default: false).
%       antisymmetry  - Logical flag for antisymmetry (default: false).
%
%   Inputs:
%       A       - A BasisFunc object.
%       options - (Optional) A structure with fields as described above.
%
%   Outputs:
%       B       - The contracted BasisFunc object after processing.
%
%   Example:
%       % Contract the rows of a BasisFunc object A with default options:
%       B = contractrow(A, struct('fast', true, 'raw', true));
%
%   See also: cleanrow, BasisFunc.introduce_coe, BasisFunc.extract_coe

    arguments
        A BasisFunc;
        options.forgetcoe = false;
        options.fast = true;
        options.hybird = false;
        options.spincoupled = false;
        options.orbcoupled = false;
        options.raw = true;
        options.sym = false;
        options.conjugate = false;
        options.antisymmetry = false;
    end

    B = A;
    if options.hybird
        error('not support yet');
        return
    end
    if iscell(A(1).BFuncL)
        error('not support yet');
        return;
    end
    A = cleanrow(A);
    if options.spincoupled
        % [Spin coupling processing can be implemented here]
    end
    if options.orbcoupled
        % [Orbital coupling processing can be implemented here]
    end
    if ~options.spincoupled && ~options.orbcoupled
        BFuncLtmp = ([A.BFuncL]);
        coeLtmp = ([A.coe]);
        BFuncLtmp = BasisFunc.introduce_coe(BFuncLtmp, coeLtmp);
        BFuncLtmp = contractrow(BFuncLtmp);
        [coeLtmp, BFuncLtmp] = BasisFunc.extract_coe(BFuncLtmp);
        switch class(BFuncLtmp)
            case 'Qnum'
                spinL = A(1).spin;
            otherwise
                spinL = A(1).spin;
        end
        orbL = A(1).BForb;
        B = BasisFunc(BFuncLtmp, spinL, 1, coeLtmp, orbL);
    end
    B = cleanrow(B);
end
