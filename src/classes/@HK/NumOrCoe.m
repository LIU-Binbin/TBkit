function [num_label, coe_label, H_hk] = NumOrCoe(H_hk)
%NUMORCOE Get or initialize numerical and coefficient labels
%
% Syntax:
%   [num_label, coe_label] = NumOrCoe(H_hk)
%   [num_label, coe_label, H_hk] = NumOrCoe(H_hk)
%
% Input:
%   H_hk - HK Hamiltonian object
%
% Outputs:
%   num_label - Numerical parameters label
%   coe_label - Coefficient parameters label
%   H_hk - Updated HK object (if requested)
%
% Description:
%   Retrieves or initializes the numerical and coefficient parameter labels:
%   1. If empty, inherits labels from TBkit superclass
%   2. If non-empty, returns stored values
%   3. Optionally updates the HK object with initialized values
%
% Note:
%   Maintains consistency between HK and TBkit parameter systems
%   Third output useful for in-place modification
%
% Example:
%   [num, coe] = Hk.NumOrCoe(); % Get current labels
%   [~, ~, Hk] = Hk.NumOrCoe(); % Initialize if empty
if isempty(H_hk.num) ||  isempty(H_hk.coe)
    [num_label,coe_label] = NumOrCoe@TBkit(H_hk);
    H_hk.num = num_label;
    H_hk.coe = coe_label;
else
    num_label = H_hk.num;
    coe_label = H_hk.coe;
end
end
