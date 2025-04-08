function fermi = fermi_get(EIGENCAR)
    % Fermi_GET - Calculate the Fermi energy from the eigenvalues of a DFT calculation
    %
    % This function computes the Fermi energy by analyzing the eigenvalues in the
    % EIGENCAR matrix. The Fermi energy is taken as the average of the
    % eigenvalues of CBM and VBM
    % corresponding to the four smallest eigenvalue differences in the middle of the list.
    %
    % Input Arguments:
    %   EIGENCAR - A matrix of eigenvalues where each column corresponds to the eigenvalues 
    %              at a particular k-point, and each row corresponds to the eigenvalue of a particular band.
    %
    % Output Arguments:
    %   fermi    - The computed Fermi energy, which is an average of the eigenvalues 
    %              of the two bands closest to the Fermi level.
    %
    % Example:
    %   fermi = fermi_get(EIGENCAR);
    %
    % Created by: [Your Name]
    
    % Get the number of Wannier functions (num_wan) and other matrix size info
    [num_wan, ~] = size(EIGENCAR);  
    
    % Sort the absolute difference between the eigenvalues of the middle two bands
    % and store the sorted indices in `label`.
    [~, label] = sort(abs(EIGENCAR(num_wan / 2 + 1, :) - EIGENCAR(num_wan / 2, :)));
     
    % Select the indices corresponding to the 4 smallest differences.
    label_list = sort(label(1:4));  
    
    % Compute the Fermi energy as the mean of the eigenvalues for the middle two bands
    % corresponding to the 3rd and 4th smallest eigenvalue differences.
    fermi = mean(mean(EIGENCAR(num_wan / 2 : num_wan / 2 + 1, label_list(3:4))));
end