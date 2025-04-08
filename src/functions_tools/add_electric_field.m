function H_E = add_electric_field(E)
% ADD_ELECTRIC_FIELD Applies an electric field to the Wannier centers and 
% returns the corresponding Hamiltonian matrix with onsite energy terms.
%
% Input:
%   E - A 1x3 vector representing the electric field in V/Å (Ex, Ey, Ez).
%       Default value is [1 0 0].
%
% Output:
%   H_E - The Hamiltonian matrix (diagonal) with onsite energy modifications
%         due to the electric field.

    % Validate input electric field vector
    arguments
        E (1,3) {mustBeReal} = [1 0 0]  % Electric field components (V/Å)
    end

    % Try to read the Wannier centers from the xyz file
    try
        % Read the data from the 'wannier90_centres.xyz' file
        xyz = readtable('wannier90_centres.xyz', 'NumHeaderLines', 2, 'FileType', 'text');
    catch
        % If file reading fails, display an error
        error('Cannot read wannier90_centres.xyz');
    end

    % Extract the Wannier centers (rows where the first column is 'X')
    wcentres = xyz(strcmp(xyz.Var1, 'X'), 2:4);
    
    % Compute the electric field interaction on the onsite energy terms
    E_onsite = table2array(wcentres) * E';  % Dot product with electric field vector
    
    % Construct the Hamiltonian matrix with onsite energy modifications
    H_E = diag(E_onsite);  % Hamiltonian is diagonal with modified onsite energies
end

% function H_E = add_electric_field(E)
% arguments
%     E (1,3) {mustBeReal} = [1 0 0] % V/Ang, [Ex Ey Ez]
% end
% try
%     xyz = readtable('wannier90_centres.xyz','NumHeaderLines',2,'FileType','text');
% catch
%     error('Can not read wannier90_centres.xyz');
% end
% wcentres = xyz( strcmp(xyz.Var1,{'X'}), 2:4);
% E_onsite = table2array(wcentres) * E';
% H_E = diag(E_onsite);
% end