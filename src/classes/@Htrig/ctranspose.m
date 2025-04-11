function H_htrig = ctranspose(H_htrig)
% CTRANSPOSE Conjugate transpose for Htrig object with dual space handling.
%   H_htrig = CTRANSPOSE(H_htrig) performs matrix transpose operations on 
%   numerical/symbolic components and handles dual space transformations 
%   for exponential-type Hamiltonians. Preserves mathematical properties 
%   in reciprocal space representations.
%
%   Inputs:
%       H_htrig - Target Htrig object to transpose
%
%   Outputs:
%       H_htrig - Modified Htrig object where:
%                 * HnumL/HcoeL: Transposed matrices for all components
%                 * Dual space transformation applied for 'exp' type
%
%   Key Operations:
%       1. Standard transpose:
%          - Numerical terms: HnumL(:,:,i) = HnumL(:,:,i)'
%          - Symbolic coefficients: HcoeL(:,:,i) = HcoeL(:,:,i)'
%       2. Dual space processing (for 'exp' type):
%          a. Invoke dualize() to generate dual representation
%          b. Reorder components using Duality_vector_dist index mapping
%
%   Example:
%       % Create exponential-type Htrig object
%       H = Htrig('Type','exp',...);
%       % Apply conjugate transpose with dual space conversion
%       H_ct = ctranspose(H);
%
%   Note:
%       - Uses non-conjugated transpose operator (') for matrix elements
%       - For 'exp' type, dualize() handles k <-> -k space transformation
%       - Original object modified due to handle class behavior
%
%   See also: dualize, transpose, ctranspose

% Determine numerical/symbolic component status
[num_label,coe_label] = H_htrig.NumOrCoe();

% Apply standard transpose to all components
for i = 1:H_htrig.Kinds
    if num_label
        H_htrig.HnumL(:,:,i) = H_htrig.HnumL(:,:,i)'; % Numerical transpose
    end
    if coe_label
        H_htrig.HcoeL(:,:,i) = H_htrig.HcoeL(:,:,i)'; % Symbolic transpose
    end
end

% Special handling for exponential-type Hamiltonians
if strcmp(H_htrig.Type,'exp')
    % Generate dual space representation
    H_htrig = H_htrig.dualize();
    
    % Reorder components using duality mapping
    if num_label
        H_htrig.HnumL(:,:,:) = H_htrig.HnumL(:,:,H_htrig.Duality_vector_dist);
    end
    if coe_label
        H_htrig.HcoeL(:,:,:) = H_htrig.HcoeL(:,:,H_htrig.Duality_vector_dist);
    end
end
end