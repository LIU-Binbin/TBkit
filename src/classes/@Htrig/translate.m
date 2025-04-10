function H_htrig = translate(H_htrig, U)
%TRANSLATE Transforms the Hamiltonian coefficients using a unitary matrix.
%
%   H_htrig = TRANSLATE(H_htrig, U) applies a unitary transformation to the Hamiltonian
%   stored in the H_htrig object using the matrix U. The transformation is applied to the
%   Hamiltonian coefficients (HcoeL and HnumL) stored within H_htrig. The function handles 
%   both double precision and symbolic unitary matrices (U).
%
%   Inputs:
%       H_htrig - An instance of the Htrig class containing the Hamiltonian and its components 
%                 (e.g., HcoeL, HnumL).
%       U       - A unitary matrix (either symbolic or numeric) used to transform the Hamiltonian.
%                 It must be square and match the dimension of the Hamiltonian components.
% 
%   Outputs:
%       H_htrig - The updated Htrig object with transformed Hamiltonian components.
%
%   Behavior:
%       - The function checks the type of the matrix U:
%         - If U is a numeric matrix (double), it applies the transformation on both HcoeL and HnumL.
%         - If U is a symbolic matrix, it only applies the transformation to HcoeL.
%       - The unitary transformation is applied as U_inv * H * U, where H is the Hamiltonian matrix.
%
%   Example:
%       % Apply a numeric transformation matrix to the Hamiltonian
%       H_htrig = translate(H_htrig, U);
%
%       % Apply a symbolic transformation matrix to the Hamiltonian
%       H_htrig = translate(H_htrig, sym('U'));
%
%   See also: Htrig, inv

    U_inv = inv(U);
    
    if isa(U,'double') % Numeric transformation
        for i = 1:H_htrig.Kinds
            H_htrig.HcoeL(:,:,i) = U_inv * H_htrig.HcoeL(:,:,i) * U;
            H_htrig.HnumL(:,:,i) = U_inv * H_htrig.HnumL(:,:,i) * U;
        end
    elseif isa(U,'sym') % Symbolic transformation
        for i = 1:H_htrig.Kinds
            H_htrig.HcoeL(:,:,i) = U_inv * H_htrig.HcoeL(:,:,i) * U;
        end
    else
        % Do nothing for unsupported types
    end
end

