function H_hk_bk = subsOper(H_hk, SymOper)
%SUBSOPER Substitute symmetry operation
%
% Syntax:
%   Hk_transformed = subsOper(Hk, SymOper)
%
% Inputs:
%   Hk - Original Hamiltonian
%   SymOper - Symmetry operation to apply
%
% Output:
%   Hk_bk - Transformed Hamiltonian
%
% Description:
%   Applies symmetry operation and returns result without modifying original
%   Essentially a non-destructive version of applyOper
%
% Note:
%   Wrapper around applyOper with different output handling
%
% Example:
%   Hk_mirror = Hk.subsOper(MirrorOperator);
H_hk_bk = H_hk.applyOper(SymOper);
end
