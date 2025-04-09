function H_hk = reseq(H_hk, basis_list)
%RESEQ Reorder basis states in Hamiltonian
%
% Syntax:
%   H_reordered = reseq(H_hk, basis_list)
%
% Inputs:
%   H_hk - Original HK Hamiltonian
%   basis_list - New basis ordering or ':' for all
%
% Output:
%   H_hk - Hamiltonian with reordered basis
%
% Description:
%   Reorders the Hamiltonian basis while preserving:
%   - Coefficient matrices (HcoeL, HnumL)
%   - Associated properties (sites, orbL, quantumL)
%   - Trigonometric terms
%
% Note:
%   Automatically initializes coefficients if empty
%   Resets Trig_to_save property
%
% Example:
%   H_new = Hk.reseq([2 4 1 3]); % Custom reordering
if isempty(H_hk.coe)||isempty(H_hk.num)
    [~,~,H_hk] = H_hk.NumOrCoe();
end
if ~isequal(basis_list,':')
    if H_hk.num
        H_hk.HnumL=H_hk.HnumL(basis_list,basis_list,:);
    else
        H_hk.HnumL = [];
    end
    if H_hk.coe
        H_hk.HcoeL=H_hk.HcoeL(basis_list,basis_list,:);
    else
        H_hk.HcoeL = sym([]);
    end
end
H_hk.Basis_num = length(basis_list);
H_hk.Trig_to_save = sym(zeros(H_hk.Basis_num));
if ~isempty(H_hk.sites)
    try
        H_hk.sites = H_hk.sites(basis_list);
    catch
    end
end
if ~isempty( H_hk.orbL )
    H_hk.orbL = H_hk.orbL(basis_list,:);
end
if ~isempty( H_hk.elementL )
    H_hk.elementL = H_hk.elementL(basis_list,:);
end
if ~isempty( H_hk.quantumL)
    H_hk.quantumL = H_hk.quantumL(basis_list,:);
end
end
