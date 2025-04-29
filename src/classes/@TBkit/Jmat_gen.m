function J_mat = Jmat_gen(H_hk,QuantumL,options)
%JMAT_GEN Generate spin matrices for Hamiltonian
%   J = Jmat_gen(H_hk, QL, opts) creates spin matrices compatible with
%   the Hamiltonian's quantum numbers.
%
%   Inputs:
%       H_hk      - HK object
%       QuantumL  - Quantum numbers (default: H_hk.quantumL)
%       options   - Optional parameters:
%           .include_0 - Include J=0 matrices (default: false)
%           .strict    - Strict compatibility check (default: false)
%
%   Output:
%       J_mat - Generated spin matrices
%
%   Example:
%       J = Jmat_gen(H, [], 'include_0', true);
arguments
H_hk  HK;
QuantumL double = H_hk.quantumL;
options.include_0 logical = false;
options.strict logical = false;
end
J_mat = Oper.spin_matrices_from_orb(QuantumL,options.include_0,options.strict);
end
