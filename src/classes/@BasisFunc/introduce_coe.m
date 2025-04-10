function [BFuncL] = introduce_coe(BFuncL, coeL)
%INTRODUCE_COE  Incorporate coefficients into each BasisFunc object.
%
%   BFuncL = INTRODUCE_COE(BFuncL, coeL) multiplies the coefficient field of each 
%   BasisFunc object in the array BFuncL by the corresponding element in the vector coeL.
%
%   Inputs:
%       BFuncL  - An array of BasisFunc objects.
%       coeL    - A vector of coefficients, with the same number of elements as BFuncL.
%
%   Output:
%       BFuncL  - The updated array of BasisFunc objects with modified 'coe' values.
%
%   Example:
%       % Given a BasisFunc array BFuncL and a coefficient vector coeL:
%       BFuncL = introduce_coe(BFuncL, coeL);
%
%   See also: extract_coe

    for i = 1:numel(BFuncL)
        BFuncL(i).coe = BFuncL(i).coe * coeL(i);
    end
end

