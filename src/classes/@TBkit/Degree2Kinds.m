function H_hk = Degree2Kinds(H_hk)
%DEGREE2KINDS Classify Hamiltonian terms by polynomial degree
%   H_hk = Degree2Kinds(H_hk) organizes Hamiltonian terms into kinds based
%   on their polynomial degree and prepares them for further processing.
%
%   Input/Output:
%       H_hk - HK object containing Hamiltonian terms
%
%   The function:
%       1. Expands polynomial terms
%       2. Classifies by degree
%       3. Normalizes and sorts terms
%       4. Performs variable substitutions
%
%   Example:
%       H = Degree2Kinds(H); % Prepare Hamiltonian for symmetry analysis
    % Extract variables up to specified dimension
    VarUsing = H_hk.VarsSeqLcart(1:H_hk.Dim);
    
    % Expand polynomial (1 + sum(VarUsing))^Degree and split terms
    str_tmp = string(expand((1 + fold(@plus, VarUsing))^H_hk.Degree));
    str_cell = strsplit(str_tmp, '+');
    
    % Store number of unique terms
    H_hk.Kinds = length(str_cell);
    
    % Initialize coefficient and degree lists
    coeff_list = zeros(H_hk.Kinds, 1);
    Degree_list = zeros(H_hk.Kinds, 1);
    
    % Extract coefficients and degrees for each term
    for i = 1:H_hk.Kinds
        tmpsym = str2sym(str_cell{i});
        Degree_list(i) = polynomialDegree(tmpsym);
        coeff_list(i) = coeffs(tmpsym);
    end
    
    % Sort terms by ascending degree
    [~, sort_label] = sort(Degree_list);
    coeff_list = coeff_list(sort_label);
    str_cell = str_cell(sort_label);
    
    % Normalize terms by coefficients and store symbolic forms
    for i = 1:H_hk.Kinds
        H_hk.HsymL(i) = str2sym(str_cell{i}) / coeff_list(i);
    end
    
    % Convert to string and perform variable substitutions
    H_hk.HstrL = string(H_hk.HsymL);
    temp_strL = H_hk.HstrL;
    temp_strL = strrep(temp_strL, 'k_x', 'x');
    temp_strL = strrep(temp_strL, 'k_y', 'y');
    temp_strL = strrep(temp_strL, 'k_z', 'z');
    
    % Apply custom substitution (assumes HK.quadk2K exists)
    for i = 1:H_hk.Kinds
        H_hk.HstrL(i) = HK.quadk2K(temp_strL(i));
    end
    
    % Store final symbolic representation
    H_hk.HsymL_xyz = str2sym(temp_strL);
end