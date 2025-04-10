function BasisFunction = rotateinner(A, abc, RightorLeft, immproper, conjugate, antisymmetry)
%ROTATEINNER  Apply inner rotation transformation to basis functions and spin.
%
%   BasisFunction = ROTATEINNER(A, abc, RightorLeft, immproper, conjugate, antisymmetry)
%   applies an inner rotation transformation to the input BasisFunc object (or a related
%   structure) A. The transformation is applied to various components of A, such as the
%   basis function list (BFuncL) and spin information. For orbitals, the transformation
%   is carried out via a dedicated method (e.g. BasisFunc.rotation_orb).
%
%   Inputs:
%       A             - A BasisFunc object or a structure containing basis function data.
%       abc           - Transformation parameter (e.g., an angle or rotation matrix) for the inner
%                       rotation.
%       RightorLeft   - A string indicating whether to apply a right or left multiplication 
%                       transformation.
%       immproper     - Logical flag indicating whether the rotation is improper (e.g., reflection).
%       conjugate     - Logical flag to indicate whether to apply conjugation during transformation.
%       antisymmetry  - Logical flag for applying antisymmetry during the rotation.
%
%   Output:
%       BasisFunction - A new BasisFunc object after applying the inner rotation to its internal 
%                       components.
%
%   Behavior:
%       - For objects of type QnumL (or similar), the function may process them directly.
%       - Otherwise, it applies rotateinner recursively to the BFuncL and spin components of A.
%       - For the orbital part, it calls a specialized rotation function (e.g., 
%         BasisFunc.rotation_orb) with appropriate parameters.
%
%   Example:
%       % Rotate the inner structure of a BasisFunc object A using given parameters:
%       BasisFunction = rotateinner(A, abc, 'right', false, true, false);
%
%   See also: BasisFunc.rotation_orb, rotateinner, QnumL
%
    switch class(A.BFuncL)
        case 'QnumL'
            % Specific processing for QnumL type can be added here.
        otherwise
            BFuncL = rotateinner(A.BFuncL, abc, RightorLeft, immproper, conjugate, antisymmetry);
            SpinL = rotateinner(A.spin, abc, RightorLeft, immproper, conjugate, antisymmetry);
            % Note: The variables R, t, and options in the following call should be defined or
            % passed as needed for the orbital rotation. Adjust as required.
            BForb = BasisFunc.rotation_orb(A.BForb, R, t, options);
            BasisFunction = BasisFunc();
            % Additional assignments to BasisFunction may be needed to combine BFuncL, SpinL, and BForb.
    end
end

