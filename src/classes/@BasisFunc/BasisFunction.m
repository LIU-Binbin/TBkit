function BasisFunction = BasisFunction(TBkitobj)
%BASISFUNCTION  Constructs a BasisFunc object (or array) from a TBkit-derived object.
%
%   BasisFunction = BASISFUNCTION(TBkitobj) creates a BasisFunc object using the
%   original basis function list, quantum numbers, and spin information extracted from
%   the input TBkit object. The input object may be of type TBkit, HR, Htrig, or HK.
%
%   Inputs:
%       TBkitobj - A TBkit-derived object representing the tight-binding system.
%
%   Outputs:
%       BasisFunction - A BasisFunc object (or an array of BasisFunc objects) that
%                       contains the basis function list, spin information, and orbital list.
%                       The reciprocal lattice matrix (Rm) from TBkitobj is assigned to each
%                       BasisFunc object.
%
%   Example:
%       % Create a BasisFunc object from a TBkit object 'tbObj':
%       BF = BasisFunction(tbObj);
%
%   See also: BasisFunc, Qnum.QnumL, Spin

    switch class(TBkitobj)
        case {'TBkit','HR','Htrig','HK'}
            [BFuncLOrigin, S, SzL] = Qnum.QnumL(TBkitobj);
            spinL = Spin(S, SzL);
            BasisFunction = BasisFunc(BFuncLOrigin, spinL, 1, 1, TBkitobj.orbL);
            for i = 1:numel(BasisFunction)
                BasisFunction(i).Rm = TBkitobj.Rm;
            end
    end
end
