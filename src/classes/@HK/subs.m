function H_hk = subs(H_hk, varargin)
%SUBS Symbolic substitution in Hamiltonian
%
% Syntax:
%   Hk = subs(Hk, old, new)
%   Hk = subs(Hk, S)       % Structure form
%   Hk = subs(Hk, old, new, ...) % Multiple pairs
%
% Inputs:
%   Hk - Original Hamiltonian
%   old - Symbolic variable/expression to replace
%   new - Replacement value
%   S - Structure with substitution pairs
%
% Output:
%   Hk - Hamiltonian with substitutions
%
% Description:
%   Performs symbolic substitution on:
%   - Coefficient matrices (HcoeL)
%   - Trigonometric terms (Trig_to_save)
%
% Note:
%   Wrapper for MATLAB's subs() with HK-specific handling
%   Supports all standard subs() syntax variants
%
% Example:
%   Hk = Hk.subs('a', 2.5) % Substitute parameter
switch length(varargin)
    case 1
        H_hk.HcoeL = subs(H_hk.HcoeL,varargin{1});
        H_hk.Trig_to_save = subs(H_hk.Trig_to_save,varargin{1});
    case 2
        H_hk.HcoeL = subs(H_hk.HcoeL,varargin{1},varargin{2});
        H_hk.Trig_to_save = subs(H_hk.Trig_to_save,varargin{1},varargin{2});
    case 3
        H_hk.HcoeL = subs(H_hk.HcoeL,varargin{1},varargin{2},varargin{3});
        H_hk.Trig_to_save = subs(H_hk.Trig_to_save,varargin{1},varargin{2},varargin{3});
end
end
