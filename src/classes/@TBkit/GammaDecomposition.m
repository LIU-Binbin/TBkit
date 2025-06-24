function [H_sym_Gamma,H_sym_Gamma_L,H_latex_Gamma] = GammaDecomposition(H_sym)
        if ~isequal(   size(H_sym) , [4,4])
                error('Gamma Decomposition requires 4*4 Ham!');
            end
            Smat_inv = gamma_matrix.S();
            Gamma_L = gamma_matrix.L();
            H_sym_L = sym(zeros(1,16));
            H_sym_Gamma_L = sym(zeros(1,16));
            %
            tmp_mat_r = real(H_sym);
            tmp_mat_i = imag(H_sym);
            H_sym_L(1:4) = diag(tmp_mat_r);
            H_sym_L(5:7) = diag(tmp_mat_r,1);
            H_sym_L(8:9) = diag(tmp_mat_r,2);
            H_sym_L(10)  = diag(tmp_mat_r,3);
            H_sym_L(11:13) = diag(tmp_mat_i,1);
            H_sym_L(14:15) = diag(tmp_mat_i,2);
            H_sym_L(16)  = diag(tmp_mat_i,3);
            %
            for i = 1:16
                Label_tmp = find(Smat_inv(i,:));
                H_sym_Gamma_L(Label_tmp) = H_sym_Gamma_L(Label_tmp)+H_sym_L(i)*Smat_inv(i,Label_tmp);
            end
            H_sym_Gamma_L = simplify(H_sym_Gamma_L);
            H_sym_Gamma = sym(0);
            H_latex_Gamma = 'H = ';
            count = 0;
            for i = 1:16
                if H_sym_Gamma_L(i)~=sym(0)
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
            %
end
