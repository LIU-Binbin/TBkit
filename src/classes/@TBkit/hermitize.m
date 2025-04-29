function H_hk = hermitize(H_hk)
%HERMITIZE Enforce Hermiticity on Hamiltonian
%   H_hk = hermitize(H_hk) ensures the Hamiltonian is Hermitian by
%   averaging with its conjugate transpose.
%
%   Input/Output:
%       H_hk - HK object to be hermitized
%
%   The operation performed is: H = (H + H')/2
%
%   Example:
%       H = hermitize(H); % Make Hamiltonian Hermitian
    % Check if numerical terms exist
    if isequal(zeros(size(H_hk.HnumL)),H_hk.HnumL)
        num_label = false;
    else
        num_label = true;
    end
    
    % Check if symbolic terms exist
    if isequal(sym(zeros(size(H_hk.HcoeL))),H_hk.HcoeL)
        coe_label = false;
    else
        coe_label = true;
    end
    
    % Hermitize components
    H_hk_bk = H_hk';
    
    if coe_label
        H_hk = (H_hk + H_hk_bk)/2;
    end
    
    if num_label
        H_hk.HnumL = (H_hk_bk.HnumL + H_hk.HnumL )/2;
    end
end