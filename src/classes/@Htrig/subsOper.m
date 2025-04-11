function H_htrig_bk = subsOper(H_htrig, SymOper)
%SUBSOPER Applies an operation to the Htrig object and returns the modified object.
%
%   H_htrig_bk = SUBSOPER(H_htrig, SymOper) applies a symmetry operation defined by the 
%   SymOper object to the input Htrig object (H_htrig) and returns a modified version 
%   of the Htrig object (H_htrig_bk).
%
%   Inputs:
%       H_htrig    - An instance of the Htrig class containing the Hamiltonian terms.
%       SymOper    - An Oper object defining the symmetry operation to be applied. 
%                   The Oper object contains information about the symmetry operation,
%                   such as the rotation matrix (R), conjugation (conjugate), antisymmetry 
%                   (antisymmetry), and any transformation matrix (U).
%
%   Outputs:
%       H_htrig_bk - An updated Htrig object after applying the symmetry operation.
%
%   Behavior:
%       - The function first checks if the Hamiltonian coefficients are all zeros; if so, 
%         it initializes the Htrig object and applies Hermitian symmetry.
%       - The function then checks the properties of the SymOper object (e.g., conjugate, 
%         antisymmetry, and rotation matrix R) to determine whether to apply the symmetry 
%         operation.
%       - If the rotation matrix is not the identity matrix, the function applies the 
%         rotation and any additional transformations (such as applying a conjugate or 
%         antisymmetry).
%       - The function returns a new Htrig object (H_htrig_bk) after the operation has been applied.
%
%   Example:
%       SymOper = Oper('R', rotMat, 'U', UMat, 'conjugate', true);
%       H_htrig_bk = subsOper(H_htrig, SymOper); % Apply symmetry operation to Htrig
%
%   See also: Oper, Htrig, applyR, applyU, hermitize

    arguments
        H_htrig Htrig;
        SymOper Oper = Oper(); % Default is an empty Oper object
    end

    if isequal(sym(zeros(size(H_htrig.HcoeL))), H_htrig.HcoeL)
        coe_label = false;
    else
        coe_label = true;
    end

    if ~coe_label
        H_htrig = H_htrig.init();
        H_htrig = H_htrig.hermitize();
    end

    if length(SymOper) == 1
        if ~SymOper.conjugate && ~SymOper.antisymmetry && isequal(SymOper.R, eye(3))
            return;
        end
        [H_htrig_bk, H_htrig] = H_htrig.applyR(inv(SymOper.R));
        if isnan(SymOper.U)
            % Additional condition when SymOper.U is NaN can be implemented here
        end
        H_htrig_bk = H_htrig_bk.applyU(SymOper.U, SymOper.conjugate, SymOper.antisymmetry);
    end
end
