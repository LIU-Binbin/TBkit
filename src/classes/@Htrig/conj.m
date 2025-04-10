function H_htrig = conj(H_htrig)
% CONJ Complex conjugate of Htrig object components.
%   H_htrig = CONJ(H_htrig) applies complex conjugation to both numerical 
%   and symbolic components of the Htrig object. Handles coefficient matrices
%   with conjugate transpose for symbolic terms to maintain Hermiticity.
%
%   Inputs:
%       H_htrig - Target Htrig object to conjugate
%
%   Outputs:
%       H_htrig - Modified Htrig object where:
%                 * HnumL(:,:,i) = conj(HnumL(:,:,i)) for numerical terms
%                 * HcoeL(:,:,i) = conj(HcoeL(:,:,i)') for symbolic coefficients
%                 * HsymL_trig = conj(HsymL_trig) for trigonometric terms
%
%   Processing Details:
%       1. Numerical components: Direct element-wise conjugation
%       2. Symbolic coefficients: Conjugate transpose to preserve matrix structure
%       3. Trigonometric symbols: Term-wise conjugation
%       4. Operates on all Kinds in the object
%
%   Example:
%       % Create Htrig object with complex elements
%       H = Htrig(...);
%       % Apply complex conjugation
%       H_conj = conj(H);
%
%   Note: 
%       - For symbolic coefficients, conjugation includes transpose (') operator
%       - Ensure proper handling of Hermitian properties in physical models
%       - Original object is modified due to handle class behavior
%
%   See also: Htrig, ctranspose, real, imag

% Get numerical/symbolic component flags
[num_label,coe_label] = H_htrig.NumOrCoe();

% Process each component type across all Kinds
for i = 1:H_htrig.Kinds
    % Numerical array conjugation
    if num_label
        H_htrig.HnumL(:,:,i) = conj(H_htrig.HnumL(:,:,i));
    end
    
    % Symbolic coefficient conjugate transpose
    if coe_label
        H_htrig.HcoeL(:,:,i) = conj(H_htrig.HcoeL(:,:,i)');
    end
end

% Conjugate trigonometric symbolic terms
H_htrig.HsymL_trig = conj(H_htrig.HsymL_trig);
end