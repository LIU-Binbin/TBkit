function [coeL, BFuncL] = extract_coe(BFuncL, options)
%EXTRACT_COE  Extract coefficients from the basis function list.
%
%   [coeL, BFuncL] = EXTRACT_COE(BFuncL) extracts the coefficient field (coe)
%   from each element of the BasisFunc object array BFuncL and returns a vector
%   coeL containing these coefficients. It also resets the coe field of each 
%   BasisFunc to one.
%
%   [coeL, BFuncL] = EXTRACT_COE(BFuncL, options) allows optional arguments:
%       options.sym      - (logical) If true, do not round the coefficients (default: false).
%       options.vpalevel - (double) Rounding precision for coefficients if not symbolic (default: 6).
%
%   Inputs:
%       BFuncL      - An array of BasisFunc objects.
%       options     - A structure with optional fields:
%                        sym      : Flag to keep symbolic coefficients (default: false).
%                        vpalevel : Precision level for rounding numeric coefficients (default: 6).
%
%   Outputs:
%       coeL        - Array of coefficients extracted from BFuncL.
%       BFuncL      - Updated BasisFunc array with the coefficient field reset to one.
%
%   Example:
%       [coefficients, newBFuncL] = extract_coe(myBasisFuncArray, struct('sym', false, 'vpalevel', 6));
%
%   See also: round

arguments
    BFuncL
    options.sym = false;
    options.vpalevel = 6;
end

coeL = ones(size(BFuncL), class(BFuncL(1).coe));
oneterm  = ones(1, 1, class(BFuncL(1).coe));
for i = 1:numel(BFuncL)
    coeL(i) = BFuncL(i).coe;
    BFuncL(i).coe = oneterm;
end
if options.sym
    % If symbolic, no rounding is performed.
else
    coeL = round(coeL, options.vpalevel);
end
end

