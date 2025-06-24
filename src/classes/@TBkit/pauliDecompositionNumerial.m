function [CoeForPauli] = pauliDecompositionNumerial(H_double)
if ~isequal(   size(H_double) , [2,2])
    error('Pauli Decomposition requires 2*2 Ham!');
end
% sigma_0+z sigma_0-z sigma_+ sigma_-
% 1/2*(sigma_0+sigma_z)
% 1/2*(sigma_0-sigma_z)
% 1/2*(sigma_x+1i*sigma_y)
% 1/2*(sigma_x-1i*sigma_y)
tmp_mat_r = real(H_double);
tmp_mat_i = imag(H_double);
S = 1/2 * [[1;0;0;1],[1;0;0;-1],[0;1;1i;0],[0;1;-1i;0]];
S = [S,S*1i];
Phi_8 = (zeros(8,1));
Phi_8(1) = tmp_mat_r(1,1);
Phi_8(2) = tmp_mat_r(2,2);
Phi_8(3) = tmp_mat_r(1,2);
Phi_8(4) = tmp_mat_r(2,1);
Phi_8(5) = tmp_mat_i(1,1);
Phi_8(6) = tmp_mat_i(2,2);
Phi_8(7) = tmp_mat_i(1,2);
Phi_8(8) = tmp_mat_i(2,1);
H_sym_Gamma_L = S * Phi_8;
count = 0;
if isa(H_double,'double')
    CoeForPauli = H_sym_Gamma_L;
else
    for i = 1:4
        if H_sym_Gamma_L(i)~=(0)
            count = count+1;
            H_sym_Gamma = H_sym_Gamma + H_sym_Gamma_L(i)*Gamma_L(i);
            str1 = latex(H_sym_Gamma_L(i));
            str2 = latex(Gamma_L(i));
            if count > 1
                H_latex_Gamma = [H_latex_Gamma,'+','\left(',str1,'\right)',str2];
            else
                H_latex_Gamma = ['\left(',str1,'\right)',str2];
            end
        end
    end
end
end
