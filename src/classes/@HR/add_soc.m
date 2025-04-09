function H_hr = add_soc(H_hr)
%ADD_SOC Add spin-orbit coupling (SOC) terms to the Hamiltonian
%   This function expands the basis to include spin degrees of freedom and 
%   adds SOC terms. It handles both numeric and symbolic coefficient formats.
%
%   Input/Output:
%       H_hr    - Input Hamiltonian without SOC. Returns Hamiltonian with SOC added.
%                 The object must have quantumL, HnumL, HcoeL, soc, coe, num properties.
%
%   Process:
%       1. If SOC is not previously enabled:
%           - Duplicate basis for spin-up/spin-down states
%           - Convert Hamiltonian to block diagonal form [H, 0; 0, H]
%       2. Update quantum numbers with spin indices (+1/-1 in 4th column)
%       3. Convert numeric coefficients to symbolic form if needed
%       4. Add atomic SOC terms using H_atom_soc() method

% Add spin dimension if SOC wasn't initialized
if ~H_hr.soc
    % Create spin-up/spin-down basis sets
    quantumList_up = H_hr.quantumL;
    quantumList_up(:,4) = 1;  % 4th column stores spin (↑=+1)
    quantumList_dn = quantumList_up;
    quantumList_dn(:,4) = -1; % Spin-down (↓=-1)
    
    % Expand Hamiltonian to block diagonal form [H, 0; 0, H]
    H_hr = blkdiag(H_hr, H_hr);  % Proper block diagonal expansion
    H_hr.soc = true;  % Mark SOC as enabled
    H_hr.quantumL = [quantumList_up; quantumList_dn]; % Combine spin bases
end

% Ensure symbolic format for SOC terms manipulation
if ~H_hr.coe
    H_hr.HcoeL = sym(H_hr.HnumL);  % Convert numeric to symbolic
    H_hr.coe = true;
    H_hr.num = false;
end

% Add atomic SOC contributions
H_hr = H_hr + H_hr.H_atom_soc(); 
end