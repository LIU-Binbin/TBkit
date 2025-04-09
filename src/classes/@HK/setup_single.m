function H_hk = setup_single(H_hk, symbolic_polynomial, i, j, silence)
%SETUP_SINGLE Configure single matrix element
%
% Syntax:
%   Hk = setup_single(Hk, poly, i, j)
%   Hk = setup_single(..., silence)
%
% Inputs:
%   Hk - HK Hamiltonian to configure
%   symbolic_polynomial - Polynomial for element (i,j)
%   i,j - Matrix indices
%   silence - Suppress output if true (default=false)
%
% Output:
%   Hk - Configured Hamiltonian
%
% Description:
%   Convenience method that:
%   1. Creates single-element matrix
%   2. Delegates to setup_rough
%   3. Simplifies single-element configuration
%
% Example:
%   Hk = Hk.setup_single(sym('a*k_x'), 1, 2);
if nargin <5
    silence = false;
end
tempmat = zeros( H_hk.Basis_num);
tempmat(i,j) =1 ;
H_hk = H_hk.setup_rough(symbolic_polynomial,tempmat,silence);
end
