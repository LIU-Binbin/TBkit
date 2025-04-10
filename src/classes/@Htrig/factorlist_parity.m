function Factorlist_parity = factorlist_parity(H_htrig)
%FACTORLIST_PARITY Calculate parity factors for Hamiltonian symmetry components
%   FACTORLIST_PARITY = FACTORLIST_PARITY(H_HTRG) computes parity transformation
%   factors for trigonometric Hamiltonian components. Returns symbolic array of
%   factors showing how components transform under k → -k inversion.
%
%INPUTS:
%   H_htrig - Htrig object containing trigonometric Hamiltonian components
%
%OUTPUTS:
%   Factorlist_parity - [N×1 sym] Array of parity factors where:
%                       H(-k) = Factorlist_parity ⊙ H(k)
%
%EXAMPLE:
%   H = Htrig(); % Initialize trigonometric Hamiltonian
%   parity_factors = factorlist_parity(H);
%
%SEE ALSO:
%   Htrig, simplify, subs

    syms k_x k_y k_z real;
    
    % Extract original Hamiltonian components
    HsymC = H_htrig.HsymL_trig;
    
    % Apply spatial inversion transformation
    HsymC = subs(HsymC, [k_x k_y k_z], -[k_x k_y k_z]);
    
    % Calculate parity factors through component-wise division
    Factorlist_parity = simplify(HsymC ./ H_htrig.HsymL_trig);
end