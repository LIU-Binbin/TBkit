function H_htrig = applyU(H_htrig,U,conjugate ,antisymmetry )
% APPLYU Applies unitary transformation and symmetry operations to Htrig components.
%   H_htrig = APPLYU(H_htrig, U, conjugate, antisymmetry) transforms the 
%   Hamiltonian coefficients/numerical terms via unitary matrix U, with optional
%   complex conjugation and antisymmetry operations. Supports direct matrix input
%   or Oper object encapsulation.
%
%   Inputs:
%       H_htrig     - Target Htrig object to transform
%       U           - Unitary matrix (n x n) or Oper object containing:
%                       * U:          Unitary matrix
%                       * conjugate:  Complex conjugation flag
%                       * antisymmetry: Antisymmetric operation flag
%       conjugate   - [Optional] Apply complex conjugation (default=false)
%       antisymmetry- [Optional] Apply antisymmetry sign flip (default=false)
%
%   Outputs:
%       H_htrig     - Transformed Htrig object with updated coefficients
%
%   Transformation Process:
%       1. Oper object handling: Extracts U/conjugate/antisymmetry if provided
%       2. Symbolic coefficient transformation:
%          a. Spatial inversion for 'sincos' type under conjugation
%          b. Apply conjugation/antisymmetry
%          c. Perform U*Coeff*U' transformation
%       3. Numerical term transformation:
%          a. Apply conjugation with parity factors
%          b. Apply antisymmetry sign
%          c. Perform batched U*Num*U' transformation
%
%   Example 1: Direct matrix input
%       U = [0 1; -1 0]; % Spin rotation
%       H_htrig = applyU(H_htrig, U, true, false);
%
%   Example 2: Oper object input
%       op = Oper('U',U,'conjugate',true);
%       H_htrig = applyU(H_htrig, op);
%
%   See also: Oper, applyR, matrixtimespage, page_mtimes_matrix

% Input validation and Oper object handling
arguments
    H_htrig Htrig;
    U = nan;
    conjugate logical = false;
    antisymmetry logical = false;
end

% Extract parameters if input is Oper object
if isa(U,'Oper')
    conjugate = U.conjugate;
    antisymmetry = U.antisymmetry;
    U = U.U;
end

% Check coefficient/numerical component status
[num_label,coe_label] = H_htrig.NumOrCoe();
U_inv = inv(U); % Precompute inverse matrix

%% Symbolic coefficient processing
if coe_label
    HcoeLtmp = H_htrig.HcoeL;
    
    % Complex conjugation handling
    if conjugate
        % Special handling for trigonometric basis
        if strcmp(H_htrig.Type,'sincos')
            % Apply spatial inversion for k -> -k transformation
            H_htrig = H_htrig.applyR(diag([-1,-1,-1]));
            HcoeLtmp = H_htrig.HcoeL;
        end
        % Apply complex conjugation
        HcoeLtmp = conj(HcoeLtmp);
    end
    
    % Antisymmetry sign flip
    if antisymmetry
        HcoeLtmp = -HcoeLtmp;
    end
    
    % Apply unitary transformation: U*Coeff*U^{-1}
    HcoeLtmp = Htrig.page_mtimes_matrix(...
        Htrig.matrix_mtimes_page(U,HcoeLtmp), U_inv);
    H_htrig.HcoeL = HcoeLtmp;
end

%% Numerical term processing
if num_label
    % Prepare batched matrices for pagemtimes
    U_page = repmat(U, [1 1 size(H_htrig.HnumL,3)]);
    U_inv_page = repmat(U_inv, [1 1 size(H_htrig.HnumL,3)]);
    
    HnumLtmp = H_htrig.HnumL;
    
    % Complex conjugation with parity factors
    if conjugate
        HnumLtmp = conj(HnumLtmp);
        HnumLtmp = Htrig.matrixtimespage(...
            H_htrig.factorlist_parity(), HnumLtmp);
    end
    
    % Antisymmetry sign flip
    if antisymmetry
        HnumLtmp = -HnumLtmp;
    end
    
    % Batched matrix multiplication: U*Num*U^{-1}
    HnumLtmp = pagemtimes(pagemtimes(U_page, HnumLtmp), U_inv_page);
    H_htrig.HnumL = HnumLtmp;
end
end