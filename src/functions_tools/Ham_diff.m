function [dH_dkx_fun, dH_dky_fun, dH_dkz_fun] = Ham_diff(Ham_obj)
% Ham_diff - Computes the derivatives of the Hamiltonian matrix with respect to k_x, k_y, and k_z.
%
% Usage:
%   [dH_dkx_fun, dH_dky_fun, dH_dkz_fun] = Ham_diff(Ham_obj)
%
% Inputs:
%   Ham_obj  - Hamiltonian object (can be of class "Htrig" or "HK").
%
% Outputs:
%   dH_dkx_fun - Function handle for the derivative of the Hamiltonian with respect to k_x.
%   dH_dky_fun - Function handle for the derivative of the Hamiltonian with respect to k_y.
%   dH_dkz_fun - Function handle for the derivative of the Hamiltonian with respect to k_z.

    % Define symbolic variables for k_x, k_y, and k_z
    syms k_x k_y k_z real;

    % Switch based on the class of the Hamiltonian object (Ham_obj)
    switch class(Ham_obj)
        case "Htrig" 
            % For Htrig class, use the corresponding Hamiltonian symmetry matrix
            symL = Ham_obj.HsymL_trig;
        case "HK"
            % For HK class, use the corresponding Hamiltonian symmetry matrix
            symL = Ham_obj.HsymL;
        otherwise
            error('Unsupported Hamiltonian class.');
    end
    
    % Compute the derivatives of the Hamiltonian with respect to k_x, k_y, and k_z
    dsym_dkx = diff(symL, k_x);  % Derivative with respect to k_x
    dsym_dky = diff(symL, k_y);  % Derivative with respect to k_y
    dsym_dkz = diff(symL, k_z);  % Derivative with respect to k_z

    % Get the number of bands (Basis_num) and the Hamiltonian matrix (HnumL)
    nbands = Ham_obj.Basis_num;
    numL = reshape(Ham_obj.HnumL, nbands^2, []);
    
    % Compute the derivatives of the Hamiltonian matrix
    dH_dkx = reshape(numL * dsym_dkx.', nbands, nbands);
    dH_dky = reshape(numL * dsym_dky.', nbands, nbands);
    dH_dkz = reshape(numL * dsym_dkz.', nbands, nbands);
    
    % Convert the symbolic derivatives into MATLAB function handles
    dH_dkx_fun = matlabFunction(dH_dkx, 'Vars', [k_x, k_y, k_z]);
    dH_dky_fun = matlabFunction(dH_dky, 'Vars', [k_x, k_y, k_z]);
    dH_dkz_fun = matlabFunction(dH_dkz, 'Vars', [k_x, k_y, k_z]);
end
