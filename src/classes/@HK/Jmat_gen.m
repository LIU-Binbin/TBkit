function J_mat = Jmat_gen(H_hk,QuantumL,options)
%JMAT_GEN Generate angular momentum matrices
%
% Syntax:
%   J_mat = Jmat_gen(H_hk)
%   J_mat = Jmat_gen(H_hk,QuantumL)
%   J_mat = Jmat_gen(H_hk,QuantumL,options)
%
% Inputs:
%   H_hk - HK object (for basis size reference)
%   QuantumL - Angular momentum quantum numbers (default=H_hk.quantumL)
%   options - Optional parameters:
%       include_0 - Include zero component if true (default=false)
%       strict - Use strict normalization if true (default=false)
%
% Output:
%   J_mat - Oper object containing Jx, Jy, Jz matrices
%
% Description:
%   Generates spin matrices for given quantum numbers using:
%   - Standard angular momentum algebra
%   - Oper.spin_matrices_from_orb() implementation
%
% See also:
%   Oper.spin_matrices_from_orb
arguments
    H_hk  HK;
    QuantumL double = H_hk.quantumL;
    options.include_0 logical = false;
    options.strict logical = false;
end
J_mat = Oper.spin_matrices_from_orb(QuantumL,options.include_0,options.strict);
end
