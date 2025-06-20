function H_hr = autohermi(H_hr, mode, options)
% AUTOHERMI Automatically enforce Hermitian symmetry on Hamiltonian
%
%   This function ensures the Hamiltonian matrix satisfies the Hermitian 
%   condition H(r) = H†(-r). It supports both symbolic ('sym') and numeric 
%   ('num') matrix modes.
%
%   Inputs:
%       H_hr    - HR object containing Hamiltonian data
%       mode    - Operation mode: 'sym' for symbolic, 'num' for numeric
%                 (default: 'sym')
%       options - Optional parameters structure:
%           enforce_list - Force conversion to list format (default: true)
%
%   Output:
%       H_hr    - Modified HR object with Hermitian-enforced Hamiltonian
%
%   Usage Example:
%       H_herm = autohermi(H_hr, 'sym', struct('enforce_list', true));
%
%   See also HR, SET_HOP, SET_HOP_MAT

    arguments
        H_hr HR;
        mode =  'sym';
        options.enforce_list = true;
    end

    % Handle input defaults
    if nargin < 2
        mode = 'sym';
    end

    % Handle format conversion
    if options.enforce_list
        if ~strcmp(H_hr.type, 'list')
            H_hr = H_hr.rewrite();
            giveback = true;
        else
            giveback = false;
        end
    else
        giveback = true;
    end

    H_hr_tmp = H_hr;  % Working copy of Hamiltonian

    switch mode
        case 'sym'  % Symbolic matrix processing
            fprintf('Strict Hermitian enforcement for symbolic Hr\n');
            fprintf('Requires real variables or explicit conj() in original Hr\n');
            
            switch H_hr.type
                case 'mat'  % Matrix format processing
                    for i = 1:H_hr.NRPTS
                        process_sym_matrix(i);
                    end
                    
                case 'list'  % List format processing
                    DIM = H_hr.Dim;
                    for i = 1:H_hr.NRPTS
                        process_sym_list(i, DIM);
                    end
                    
                otherwise
                    error('Unsupported format for sym mode');
            end

        case 'num'  % Numeric matrix processing
            fprintf('Relaxed Hermitian enforcement for numeric Hr\n');
            for i = 1:H_hr.NRPTS
                process_num_matrix(i);
            end

        otherwise
            error('Invalid mode. Use ''sym'' or ''num''');
    end

    % Final conversion if needed
    if giveback
        H_hr = H_hr_tmp.rewind();
    else
        H_hr = H_hr_tmp;
    end

    % Nested processing functions
    function process_sym_matrix(i)
        % Process symbolic matrix format
        vector_tmp = H_hr.vectorL(i,:);
        vector_tmp_oppo = -vector_tmp;
        [~, j] = ismember(vector_tmp_oppo, H_hr_tmp.vectorL, 'rows');
        
        fprintf('Checking NRPT %d:\n Vector: [%d %d %d]\n Opposite: [%d %d %d] at NRPT %d\n',...
                i, vector_tmp, vector_tmp_oppo, j);
        
        if i == j  % Diagonal term (home cell)
            if ~isequal(H_hr.HcoeL(:,:,i), H_hr_tmp.HcoeL(:,:,j)')
                enforce_homecell_hermi(i);
            end
        else
            handle_offdiagonal_term(i, j, vector_tmp, vector_tmp_oppo);
        end
    end

    function process_sym_list(i, DIM)
        % Process symbolic list format
        vector_tmp = H_hr.vectorL(i,:);
        vector_tmp_oppo = zeros(1, DIM+2);
        vector_tmp_oppo(1:DIM) = -vector_tmp(1:DIM);
        vector_tmp_oppo(DIM+1) = vector_tmp(DIM+2);
        vector_tmp_oppo(DIM+2) = vector_tmp(DIM+1);
        [~, j] = ismember(vector_tmp_oppo, H_hr_tmp.vectorL, 'rows');
        
        if i == j  % Diagonal term
            if ~isequal(H_hr.HcoeL(i), H_hr_tmp.HcoeL(j)')
                handle_homecell_nonhermi(i, j, vector_tmp, DIM);
            end
        else
            handle_list_offdiagonal(i, j, vector_tmp, vector_tmp_oppo, DIM);
        end
    end

    function process_num_matrix(i)
        % Process numeric matrix format
        vector_tmp = H_hr.vectorL(i,:);
        vector_tmp_oppo = -vector_tmp;
        [~, j] = ismember(vector_tmp_oppo, H_hr_tmp.vectorL, 'rows');
        
        fprintf('Checking NRPT %d:\n Vector: [%d %d %d]\n Opposite: [%d %d %d] at NRPT %d\n',...
                i, vector_tmp, vector_tmp_oppo, j);
        
        if i == j  % Diagonal term
            if ~isequal(H_hr.HnumL(:,:,i), H_hr_tmp.HnumL(:,:,j)')
                enforce_numeric_hermi(i, j);
            end
        else
            handle_offdiagonal_num(i, j, vector_tmp);
        end
    end

    % Helper functions
    function enforce_homecell_hermi(i)
        fprintf('Non-Hermitian homecell Hamiltonian found at NRPT %d\n', i);
        disp('Original matrix:');
        disp(H_hr.HcoeL(:,:,i));
        corrected = (H_hr.HcoeL(:,:,i)' + H_hr_tmp.HcoeL(:,:,i))/2;
        H_hr_tmp = H_hr_tmp.set_hop_mat(corrected, [0,0,0], 'sym');
        disp('Corrected matrix:');
        disp(H_hr_tmp.HcoeL(:,:,i));
    end

    function handle_offdiagonal_term(i, j, vec, vec_oppo)
        if j == 0
            fprintf('Creating missing opposite vector at NRPT %d\n', i);
            H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr.HcoeL(:,:,i)', vec_oppo, 'sym');
        else
            compare_and_correct(i, j, vec, vec_oppo);
        end
    end

    function compare_and_correct(i, j, vec, vec_oppo)
        if ~isequal(H_hr.HcoeL(:,:,i), H_hr_tmp.HcoeL(:,:,j)')
            fprintf('Hermitian mismatch between NRPT %d and %d\n', i, j);
            n1 = nnz(H_hr.HcoeL(:,:,i));
            n2 = nnz(H_hr_tmp.HcoeL(:,:,j)');
            
            if n1 >= n2
                fprintf('Enforcing stronger term from NRPT %d\n', i);
                H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr.HcoeL(:,:,i)', vec_oppo, 'sym');
            else
                fprintf('Enforcing stronger term from NRPT %d\n', j);
                H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr_tmp.HcoeL(:,:,j)', vec, 'sym');
            end
        else
            disp('Hermitian condition satisfied');
        end
    end

    % Additional nested helper functions for list processing
    function handle_homecell_nonhermi(i, j, vec, DIM)
        if H_hr.HcoeL(i) == sym(0)
            fprintf('Adding missing homecell term from NRPT %d\n', j);
            H_hr_tmp = H_hr_tmp.set_hop(H_hr.HcoeL(j)', vec(DIM+1), vec(DIM+2), vec(1:DIM), 'sym');
        elseif H_hr.HcoeL(j) == sym(0)
            fprintf('Adding missing homecell term from NRPT %d\n', i);
            H_hr_tmp = H_hr_tmp.set_hop(H_hr.HcoeL(i)', vec(DIM+1), vec(DIM+2), vec(1:DIM), 'sym');
        else
            fprintf('Averaging non-Hermitian homecell terms\n');
            tmpsym = (H_hr.HcoeL(i) + H_hr.HcoeL(j)')/2;
            H_hr_tmp = H_hr_tmp.set_hop(tmpsym, vec(DIM+1), vec(DIM+2), vec(1:DIM), 'sym');
        end
    end

    function handle_list_offdiagonal(i, j, vec, vec_oppo, DIM)
        if j == 0
            fprintf('Creating missing opposite vector for %s \n', num2str(vec));
            H_hr_tmp = H_hr_tmp.set_hop(H_hr.HcoeL(i)', vec_oppo(DIM+1), vec_oppo(DIM+2), vec_oppo(1:DIM), 'sym');
        else
            if ~isequal(H_hr.HcoeL(i), H_hr_tmp.HcoeL(j)')
                fprintf('Resolving Hermitian mismatch between NRPT %s and %s\n', num2str(vec), num2str(vec_oppo));
                tmpsym = (H_hr.HcoeL(i) + H_hr.HcoeL(j)')/2;
                H_hr_tmp = H_hr_tmp.set_hop(tmpsym, vec(DIM+1), vec(DIM+2), vec(1:DIM), 'sym');
                H_hr_tmp = H_hr_tmp.set_hop(tmpsym', vec_oppo(DIM+1), vec_oppo(DIM+2), vec_oppo(1:DIM), 'sym');
            end
        end
    end

    function enforce_numeric_hermi(i, j)
        fprintf('Enforcing Hermiticity on numeric homecell at NRPT %d\n', i);
        corrected = (H_hr.HnumL(:,:,i)' + H_hr_tmp.HnumL(:,:,j))/2;
        H_hr_tmp = H_hr_tmp.set_hop_mat(corrected, [0,0,0], 'set');
    end

    function handle_offdiagonal_num(i, j, vec)
        if j == 0
            fprintf('Creating missing numeric opposite vector at NRPT %d\n', i);
            H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr.HnumL(:,:,i)', -vec, 'set');
        else
            if ~isequal(H_hr.HnumL(:,:,i), H_hr_tmp.HnumL(:,:,j)')
                fprintf('Correcting numeric Hermitian mismatch between NRPT %d and %d\n', i, j);
                n1 = nnz(H_hr.HnumL(:,:,i));
                n2 = nnz(H_hr_tmp.HnumL(:,:,j)');
                if n1 >= n2
                    H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr.HnumL(:,:,i)', -vec, 'set');
                else
                    H_hr_tmp = H_hr_tmp.set_hop_mat(H_hr_tmp.HnumL(:,:,j)', vec, 'set');
                end
            end
        end
    end

end