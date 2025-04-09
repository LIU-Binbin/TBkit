function [sym_coe,sym_single,str_single] = coeff_extract(sym_term,Dim)
%COEFF_EXTRACT Extract coefficients and terms from symbolic expression
%
% Syntax:
%   [sym_coe,sym_single] = coeff_extract(sym_term)
%   [sym_coe,sym_single,str_single] = coeff_extract(sym_term,Dim)
%
% Inputs:
%   sym_term - Symbolic expression to analyze
%   Dim - Dimension of k-space (default=3)
%
% Outputs:
%   sym_coe - Coefficients of each term
%   sym_single - Symbolic terms
%   str_single - String representations of terms
%
% Description:
%   Decomposes a symbolic expression into its constituent terms and
%   coefficients with respect to momentum variables k_x, k_y, k_z.
%   Useful for Hamiltonian term analysis and manipulation.
%
% Example:
%   [c,t,s] = coeff_extract(sym('a*k_x + b*k_y^2'),2)
arguments
    sym_term
    Dim = 3;
end
VarsSeqLcart = [sym('k_x'),sym('k_y'),sym('k_z'),sym('k_w')];
[sym_coe,sym_single] = coeffs(sym_term,VarsSeqLcart(1:Dim));
str_single = string(sym_single);
end
