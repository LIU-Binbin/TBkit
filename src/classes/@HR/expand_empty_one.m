function H_hr = expand_empty_one(H_hr,orbOne,QuantumOne,elementOne)
% EXPAND_EMPTY_ONE Add empty orbitals to HR object
%
%   H_hr = EXPAND_EMPTY_ONE(H_hr,orbOne,QuantumOne,elementOne) adds empty
%   orbitals to the Hamiltonian with specified properties.
%
%   INPUT ARGUMENTS:
%       H_hr - Hamiltonian in HR format
%       orbOne - Orbital positions (default: zeros(1,Dim))
%       QuantumOne - Quantum numbers (default: [1,0,0,1])
%       elementOne - Element identifiers (default: 1)
%
%   OUTPUT ARGUMENTS:
%       H_hr - Modified Hamiltonian with added empty orbitals
%
%   NOTES:
%       - Expands both orbital list and Hamiltonian matrices
%       - Handles both 'list' and 'mat' Hamiltonian types
%       - Initializes new elements to zero
%
%   SEE ALSO:
%       HR
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

arguments
H_hr;
orbOne = zeros(1,H_hr.Dim);
QuantumOne = [1,0,0,1];
elementOne = 1;
end
H_hr.orbL = [H_hr.orbL;orbOne];
H_hr.quantumL = [H_hr.quantumL;QuantumOne];
H_hr.elementL = [H_hr.elementL;elementOne];
if strcmp(H_hr.Type ,'list')
NRPTS_new = H_hr.NRPTS +size(orbOne,1);
H_hr.HcoeL(NRPTS_new,1) = sym(0);
H_hr.HnumL(NRPTS_new,1) = 0;
else
WANNUM = H_hr.WAN_NUM+size(orbOne,1);
H_hr.HcoeL(WANNUM,WANNUM,:) = sym(0);
H_hr.HnumL(WANNUM,WANNUM,:) = 0;
end
end
