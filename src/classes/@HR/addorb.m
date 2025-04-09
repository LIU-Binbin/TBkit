function H_hr = addorb(H_hr, orblist, options)
%ADDORB Expand Hamiltonian by adding duplicated orbitals
%   This function extends the basis space by replicating specified orbitals,
%   creating block-structured Hamiltonian couplings. Handles numeric coefficient format.
%
%   Inputs:
%       H_hr    - Original HR object containing Hamiltonian
%       orblist - List of orbital indices to duplicate (1-based indices)
%       options.inner - [Optional] Flag for intra-block coupling (Reserved for future use)
%
%   Output:
%       H_hr    - Modified HR object with expanded Hamiltonian matrix
%
%   Features:
%       - Creates (N+norb)x(N+norb) block matrix from original NxN matrix
%       - Preserves original orbital ordering, appends new orbitals at end
%       - Copies row/column elements from specified orbitals for coupling blocks
%       - Maintains third dimension (e.g., R-vectors) unchanged

% Argument validation with type constraints
arguments
    H_hr HR;
    orblist double {mustBePositive, mustBeInteger};
    options.inner logical = true;
end

% Cache original Wannier function count
orig_wan_num = H_hr.WAN_NUM;
norb = length(orblist);

% Create expanded Hamiltonian container
HnumLtmp = H_hr.HnumL;  % Original Hamiltonian data
new_size = size(HnumLtmp) + [norb, norb, 0];
HnumLtmp2 = zeros(new_size);

% Preserve original block in top-left quadrant
HnumLtmp2(1:orig_wan_num, 1:orig_wan_num, :) = HnumLtmp;

% Add orbital coupling blocks:
% Bottom-left block: new orbitals × original basis
HnumLtmp2(orig_wan_num+1:end, 1:orig_wan_num, :) = HnumLtmp(orblist, :, :);

% Top-right block: original basis × new orbitals
HnumLtmp2(1:orig_wan_num, orig_wan_num+1:end, :) = HnumLtmp(:, orblist, :);

% Update Hamiltonian storage
H_hr.HnumL = HnumLtmp2;

% Note: Should update H_hr.WAN_NUM += norb and handle quantumL in actual implementation
% This prototype only demonstrates matrix expansion logic
end