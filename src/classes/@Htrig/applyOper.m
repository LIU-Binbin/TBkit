function H_htrig = applyOper(H_htrig,SymOper,options)
% APPLYOPER Applies symmetry operator to Htrig object and enforces constraints.
%   H_htrig = APPLYOPER(H_htrig, SymOper, options) transforms the Hamiltonian
%   according to the specified symmetry operator and resolves resulting constraints.
%
%   Inputs:
%       H_htrig    - Target Htrig object to operate on
%       SymOper    - Oper object(s) defining symmetry operation(s). Contains:
%                      * R:        3x3 rotation matrix
%                      * U:        Unitary transformation matrix
%                      * conjugate: Complex conjugation flag
%                      * antisymmetry: Antisymmetric operation flag
%       options.fast - [Optional] Fast mode skips intermediate checks (default=true)
%
%   Outputs:
%       H_htrig    - Transformed Htrig object satisfying symmetry constraints
%
%   Operation Flow:
%       1. Coefficient initialization: Ensures coefficient matrices exist
%       2. Single operator handling:
%          a. Apply inverse rotation via applyR()
%          b. Apply unitary/conjugation via applyU()
%          c. Derive real/imaginary constraints from symmetry
%          d. Substitute solutions into coefficient matrices
%       3. Multi-operator handling: Recursively apply individual operators
%
%   Example:
%       % Create rotation operator and apply to Htrig
%       R = rotx(90); % 90-degree rotation about x-axis
%       SymOp = Oper('R',R,'U',eye(2));
%       H_htrig = applyOper(H_htrig, SymOp);
%
%   See also: Oper, applyR, applyU, hermitize

% Input validation and type checking
arguments
    H_htrig Htrig;
    SymOper Oper = Oper();
    options.fast = true;
end

% Ensure coefficient matrices are initialized
[~,coe_label] = H_htrig.NumOrCoe();
if ~coe_label
    H_htrig = H_htrig.init();
    H_htrig = H_htrig.hermitize();
end

% Main operator processing
if length(SymOper) == 1  % Single operator case
    % Skip trivial operations
    if ~SymOper.conjugate && ~SymOper.antisymmetry && isequal(SymOper.R,eye(3))
        return;
    end
    
    % Apply inverse rotation transformation
    [H_htrig_bk,H_htrig]  = H_htrig.applyR(inv(SymOper.R));
    
    % Apply unitary/antisymmetry/conjugation operations
    H_htrig_bk  = H_htrig_bk.applyU(SymOper.U,SymOper.conjugate,SymOper.antisymmetry);
    
    % Derive symmetry constraint equations
    Equationlist_r = (real(H_htrig.HcoeL - H_htrig_bk.HcoeL) == 0);
    Equationlist_i = (imag(H_htrig.HcoeL - H_htrig_bk.HcoeL) == 0);
    
    % Solve and substitute constraints
    Equationlist_r = Htrig.isolateAll(Equationlist_r);
    Equationlist_i = Htrig.isolateAll(Equationlist_i);
    HcoeLtmp = H_htrig.HcoeL ;
    HcoeLtmp_r = subs(real(HcoeLtmp),lhs(Equationlist_r),rhs(Equationlist_r));
    HcoeLtmp_i = subs(imag(HcoeLtmp),lhs(Equationlist_i),rhs(Equationlist_i));
    H_htrig.HcoeL = HcoeLtmp_r + 1i*HcoeLtmp_i;
else  % Multiple operators - recursive processing
    for i = 1:length(SymOper)
        H_htrig = H_htrig.applyOper(SymOper(i));
    end
end
end