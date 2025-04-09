function H_hr = charalize(H_hr)
% CHARALIZE Zero off-diagonal couplings between same-element orbitals
%
%   H_hr = CHARALIZE(H_hr) eliminates intra-element orbital couplings by
%   setting all off-diagonal Hamiltonian terms between orbitals of the
%   same atomic element to zero. This operation is typically used to
%   enforce specific symmetry constraints or simplify tight-binding models.
%
%   Input:
%       H_hr - HR object containing Hamiltonian data with properties:
%           elementL: Atom element labels (cellstr or char array)
%           HnumL:    Numeric Hamiltonian matrices (WAN_NUM x WAN_NUM x NRPTS)
%           WAN_NUM:  Number of Wannier orbitals
%   Output:
%       H_hr - Modified HR object with zeroed intra-element couplings
%
%   Example:
%       H_hr = charalize(H_hr); % Remove same-element orbital couplings
%
%   See also HR, AUTOHERMI

    % Determine element comparison method based on data type
    if iscellstr(H_hr.elementL)
        compare = @(i,j) strcmp(H_hr.elementL{i}, H_hr.elementL{j});
    elseif ischar(H_hr.elementL)
        elementCell = cellstr(H_hr.elementL);
        compare = @(i,j) strcmp(elementCell{i}, elementCell{j});
    else
        error('Element labels must be cell array of strings or character array');
    end

    % Zero out off-diagonal terms for same-element orbitals
    for i = 1:H_hr.WAN_NUM
        for j = 1:H_hr.WAN_NUM
            if compare(i,j) && i ~= j
                H_hr.HnumL(i,j,:) = 0; % Operates on all k-points/R-vectors
            end
        end
    end
end