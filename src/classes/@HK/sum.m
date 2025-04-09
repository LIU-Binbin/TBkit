function H_hk = sum(H_hk_list)
%SUM Sum array of HK Hamiltonians
%
% Syntax:
%   Hk_total = sum(Hk_array)
%
% Input:
%   Hk_list - Array of HK objects
%
% Output:
%   Hk - Combined Hamiltonian
%
% Description:
%   Iteratively sums Hamiltonians using plus() operator
%   Preserves all properties of first Hamiltonian
%
% Note:
%   Requires consistent basis_num and Degree
%   Uses overloaded + operator internally
%
% Example:
%   Hk_total = sum([Hk1, Hk2, Hk3]);
H_hk = H_hk_list(1);
for i = 2:length(H_hk_list)
    H_hk = H_hk + H_hk_list(i);
end
end
