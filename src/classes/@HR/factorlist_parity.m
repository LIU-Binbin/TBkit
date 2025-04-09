function Factorlist_parity = factorlist_parity(H_hr)
% FACTORLIST_PARITY Calculate parity factors for Hamiltonian
%
%   Factorlist_parity = FACTORLIST_PARITY(H_hr) computes parity factors
%   by comparing Hamiltonian at k and -k points.
%
%   INPUT ARGUMENTS:
%       H_hr - Hamiltonian in HR format with symbolic HsymL_trig
%
%   OUTPUT ARGUMENTS:
%       Factorlist_parity - Symbolic expression of parity factors
%
%   NOTES:
%       - Uses symbolic variables k_x, k_y, k_z
%       - Returns simplified ratio H(-k)/H(k)
%
%   SEE ALSO:
%       HR
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

syms k_x k_y k_z real;
HsymC = H_hr.HsymL_trig;
HsymC = subs(HsymC,[k_x k_y k_z],-[k_x k_y k_z]);
Factorlist_parity = simplify(HsymC./H_hr.HsymL_trig);
end
