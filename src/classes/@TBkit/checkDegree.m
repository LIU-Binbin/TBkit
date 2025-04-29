function Degree = checkDegree(Hsym,Dim)
%CHECKDEGREE Calculate maximum polynomial degree of Hamiltonian terms
%   Degree = checkDegree(Hsym, Dim) determines the maximum polynomial degree
%   of symbolic Hamiltonian terms in specified dimensions.
%
%   Inputs:
%       Hsym - Symbolic Hamiltonian expression
%       Dim  - Number of dimensions to consider (default: 3)
%
%   Output:
%       Degree - Maximum polynomial degree found
%
%   Example:
%       max_degree = checkDegree(H, 2); % Check in 2D
arguments
Hsym
Dim = 3;
end
VarsSeqLcart = [sym('k_x'),sym('k_y'),sym('k_z'),sym('k_w')];
Degree_mat = polynomialDegree(Hsym,VarsSeqLcart(1:Dim));
Degree = max(Degree_mat,[],'all');
end
