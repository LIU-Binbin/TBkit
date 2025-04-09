function Degree = checkDegree(Hsym,Dim)
%CHECKDEGREE Determine polynomial degree of symbolic Hamiltonian
%
% Syntax:
%   Degree = checkDegree(Hsym)
%   Degree = checkDegree(Hsym,Dim)
%
% Inputs:
%   Hsym - Symbolic Hamiltonian matrix
%   Dim - Dimension of k-space (default=3)
%
% Output:
%   Degree - Highest polynomial degree found in Hamiltonian
%
% Description:
%   Calculates the maximum polynomial degree of a symbolic Hamiltonian with
%   respect to momentum variables k_x, k_y, k_z (and optionally k_w).
%   Used to determine the order of k·p expansion.
%
% Example:
%   Degree = checkDegree(sym('k_x^2 + k_y*k_z'),3) % Returns 2
arguments
    Hsym
    Dim = 3;
end
VarsSeqLcart = [sym('k_x'),sym('k_y'),sym('k_z'),sym('k_w')];
Degree_mat = polynomialDegree(Hsym,VarsSeqLcart(1:Dim));
Degree = max(Degree_mat,[],'all');
end
