function C = eq(A, B, options)
%EQ  Overloaded equality operator for BasisFunc objects.
%
%   C = EQ(A, B) compares two BasisFunc objects A and B for equality.
%
%   C = EQ(A, B, options) performs the equality check with additional options.
%
%   Inputs:
%       A           - A BasisFunc object.
%       B           - A BasisFunc object.
%       options     - A structure with the following fields:
%                       strict   - (logical) If true, requires strict equality (default: false).
%                       spin     - (logical) If true, compares the spin properties (default: false).
%                       orb      - (logical) If true, compares the orbital information (BForb) (default: true).
%                       BFuncL   - (logical) If true, compares the basis function list (BFuncL) (default: true).
%
%   Output:
%       C           - Boolean value; true if A and B are considered equal based on the criteria, 
%                     false otherwise.
%
%   Behavior:
%       - The function first checks whether the class of the BFuncL property of A and B are equal.
%       - If the 'spin' option is enabled, it further checks the equality of the spin properties.
%       - If the 'orb' option is enabled, it compares the orbital information (BForb) using Oper.allclose.
%       - Finally, if the 'BFuncL' option is enabled, it compares the basis function lists.
%
%   Example:
%       % Compare two BasisFunc objects with orbital and basis function list checks:
%       C = eq(BasisFunc1, BasisFunc2, struct('strict', true, 'spin', false, 'orb', true, 'BFuncL', true));
%
%   See also: BasisFunc, Oper.allclose

arguments
    A BasisFunc;
    B BasisFunc;
    options.strict = false;
    options.spin = false;
    options.orb = true;
    options.BFuncL = true;
end

C = true;
if class(A.BFuncL) ~= class(B.BFuncL)
    C = false;
    return;
end
if options.spin
    if ~eq(A.spin, B.spin, 'strict', options.strict)
        C = false;
        return;
    end
end
if options.orb
    if isempty(A.BForb) && ~isempty(B.BForb)
        C = false;
        return;
    end
    if ~isempty(A.BForb) && isempty(B.BForb)
        C = false;
        return;
    end
    if ~Oper.allclose(A.BForb, B.BForb)
        C = false;
        return;
    end
end
if options.BFuncL
    if ~eq(A.BFuncL, B.BFuncL, 'strict', options.strict)
        C = false;
        return;
    end
end
end
