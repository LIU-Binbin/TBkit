function C = innertimes(A, B, options)
%INNERTIMES  Compute the inner (element-wise) product of two BasisFunc objects.
%
%   C = INNERTIMES(A, B) computes the inner product between two BasisFunc objects A and B.
%
%   C = INNERTIMES(A, B, options) computes the inner product with additional options:
%       options.forgetcoe   - (logical) If true, the coefficients may be ignored (default: false).
%       options.fast        - (logical) If true, uses a fast computation mode (default: true).
%       options.hybird      - (logical) Hybrid mode flag; if true, uses hybrid contraction (default: false).
%       options.spincoupled - (logical) If true, considers spin coupling in the computation (default: false).
%       options.orbcoupled  - (logical) If true, considers orbital coupling in the computation (default: false).
%       options.raw         - (logical) If true, returns a raw result without further processing (default: true).
%
%   Inputs:
%       A         - A BasisFunc object.
%       B         - A BasisFunc object.
%       options   - (Optional) A structure with the options described above.
%
%   Output:
%       C         - A new BasisFunc object representing the inner product of A and B.
%
%   Example:
%       C = innertimes(A, B, struct('fast', true, 'raw', true));
%
%   See also: extract_coe, BasisFunc

    arguments
        A BasisFunc;
        B BasisFunc;
        options.forgetcoe = false;
        options.fast = true;
        options.hybird = false;
        options.spincoupled = false;
        options.orbcoupled = false;
        options.raw = true;
    end

    optionsCell = namedargs2cell(options);
    if ~options.hybird && ~options.spincoupled && ~options.orbcoupled
        coeLtmp0 = A.coe * B.coe;
        BFuncLtmp = (A.BFuncL .* B.BFuncL);
        [coeLtmp, BFuncLtmp] = BasisFunc.extract_coe(BFuncLtmp);
        spinL = A(1).spin;
        orbL = A(1).BForb;
        C = BasisFunc(BFuncLtmp, spinL, 1, coeLtmp0 * coeLtmp, 'orbL', orbL);
    end
end

