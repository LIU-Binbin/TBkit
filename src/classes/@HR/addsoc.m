function H_hr = addsoc(H_hr, quantumL)
%ADDSOC Add atomic spin-orbit coupling (SOC) to the Hamiltonian
%   This function incorporates atomic-level SOC terms into the Hamiltonian.
%   Allows optional update of quantum number descriptors for SOC parameters.
%
%   Inputs:
%       H_hr     - HR object containing the base Hamiltonian
%       quantumL - (Optional) Updated quantum number list [Nx4 matrix]
%                  Columns: [n, l, j, m_j] atomic orbital descriptors
%
%   Output:
%       H_hr     - Modified HR object with SOC terms added
%
%   Features:
%       - Preserves existing Hamiltonian terms
%       - Atomic SOC terms are added through H_atom_soc() method
%       - Updates quantum numbers if provided
%       - Maintains original coefficient type (numeric/symbolic)

% Handle optional quantum number update
if nargin > 1
    % Validate quantum number matrix dimensions
    assert(size(quantumL,2) >= 4, 'quantumL must have at least 4 columns');
    H_hr.quantumL = quantumL;
end

% Add atomic SOC terms using class method
H_hr = H_hr + H_hr.H_atom_soc();

% Auto-enable SOC flag if exists in object
if isprop(H_hr, 'soc')
    H_hr.soc = true;
end
end