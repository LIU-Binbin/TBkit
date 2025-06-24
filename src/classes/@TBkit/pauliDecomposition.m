     function [H_sym_pauli,H_sym_pauli_L,H_latex_pauli] = pauliDecomposition(H_sym)
            if isa(H_sym,'vasplib') || isa(H_sym,'HR') || isa(H_sym,'Htrig') || isa(H_sym,'HK')
                H_sym = H_sym.sym();
            elseif isa(H_sym,'double')
                H_sym = sym(H_sym);
            elseif isa(H_sym,'sym')
            else
                try
                    H_sym = sym(H_sym);
                catch ME
                    error('cant convert into symbolic');
                end
            end
            if ~isequal(   size(H_sym) , [4,4]) && ~isequal(   size(H_sym) , [2,2])
                error('Pauli Decomposition requires 4*4/2*2 Ham!');
            end
            if isequal(   size(H_sym) , [4,4])
                [~,H_sym_pauli_L] = TBkit.GammaDecomposition(H_sym);
                Pauli_L = gamma_matrix.pauli_L();
                
                H_sym_pauli = sym(0);
                H_latex_pauli = 'H = ';
                count = 0;
                ChooseL = logical(1:16);
                for i = 1:16
                    if H_sym_pauli_L(i)~=sym(0)
                        count = count+1;
                        H_sym_pauli = H_sym_pauli + H_sym_pauli_L(i)*Pauli_L(i);
                        str1 = latex(H_sym_pauli_L(i));
                        str2 = latex(Pauli_L(i));
                        if count >1
                            H_latex_pauli = [H_latex_pauli,'+','\left(',str1,'\right)',str2];
                        else
                            H_latex_pauli = ['\left(',str1,'\right)',str2];
                        end
                    else
                        ChooseL(i) = false;
                    end
                end
                H_sym_pauli_L = [H_sym_pauli_L;Pauli_L];
                H_sym_pauli_L = H_sym_pauli_L(:,ChooseL);
            end
            if isequal(   size(H_sym) , [2,2])
                tmp_mat_r = real(H_sym);
                tmp_mat_i = imag(H_sym);
                H_sym_L = sym(zeros(1,4));
                H_sym_pauli_L = sym(zeros(1,4));
%                 sigma_L = pauli_matrix();
                %
                Smat = [
                    1 0 0 1;
                    1 0 0 -1;
                    0 1 1i 0;
                    0 1 -1i 0];
%                 Smat = sym([1/2,1/2,0,0; ...
%                                 0,0,1,0; ...
%                                 0,0,0,-1; ...
%                                 1/2,-1/2,0,0]);
                Smat_inv = inv(Smat);
                %
                syms sigma_0 sigma_x sigma_y sigma_z real;
                Pauli_L = [sigma_0 sigma_x sigma_y sigma_z];
                H_sym_L(1:2) = diag(H_sym);
                H_sym_L(3) = diag(H_sym,1);
                H_sym_L(4) = diag(H_sym,-1);
                for i = 1:numel(H_sym_L)
                    Label_tmp = find(Smat_inv(i,:));
                    H_sym_pauli_L(i) = H_sym_L(Label_tmp)*Smat_inv(i,Label_tmp)';
                end
%                 H_sym_pauli_L_temp = Smat*H_sym_L';
%                 H_sym_pauli_L = H_sym_pauli_L_temp';
                H_sym_pauli_L = simplify(H_sym_pauli_L);
                H_sym_pauli = sym(0);
                H_latex_pauli = 'H = ';
                count = 0;
                ChooseL = logical(1:4);
                %Pauli_L = sym(Pauli_L);
                for i = 1:4
                    if H_sym_pauli_L(i)~=sym(0)
                        count = count+1;
                        H_sym_pauli = H_sym_pauli + H_sym_pauli_L(i)*sym(Pauli_L(i));
                        str1 = latex(H_sym_pauli_L(i));
                        str2 = latex(sym(Pauli_L(i)));
                        if count >1
                            H_latex_pauli = [H_latex_pauli,'+','\left(',str1,'\right)',str2];
                        else
                            H_latex_pauli = ['\left(',str1,'\right)',str2];
                        end
                    else
                        ChooseL(i) = false;
                    end
                end
                H_sym_pauli_L = [H_sym_pauli_L;Pauli_L];
                H_sym_pauli_L = H_sym_pauli_L(:,ChooseL);
            end
        end
   