function H_hr = hr_converter(mode, write)
% hr_converter - Converts the Hamiltonian matrix between old and new formats.
%
% Usage:
%   H_hr = hr_converter(mode, write)
%
% Inputs:
%   mode   - Specifies the conversion mode. Options:
%            'old_to_new' : Converts from old to new format.
%            'new_to_old' : Converts from new to old format.
%   write  - Logical flag. If true, the converted Hamiltonian is saved to a file.
%
% Outputs:
%   H_hr   - The Hamiltonian object with the converted matrix.
%
% Example:
%   H_hr = hr_converter('old_to_new', true);

    % Validate input arguments using MATLAB's argument validation
    arguments
        mode {mustBeMember(mode, {'old_to_new', 'new_to_old'})} = 'old_to_new';
        write logical = true;
    end

    % Load the Hamiltonian object from the 'wannier90' file
    H_hr = HR.from_wannier90();
    
    % Initialize a temporary matrix to store the converted Hamiltonian
    HnumL_tmp = zeros(size(H_hr.HnumL));
    
    % Calculate the number of orbitals (half the number of bands)
    Norb = H_hr.Nbands / 2;

    % Perform conversion based on the specified mode
    switch mode
        case 'old_to_new'
            % Convert from 'old' format to 'new' format
            for i = 1:Norb
                for j = 1:Norb
                    % Mapping elements from old to new format
                    HnumL_tmp(2*i-1, 2*j-1, :) = H_hr.HnumL(i, j, :);
                    HnumL_tmp(2*i, 2*j-1, :) = H_hr.HnumL(i + Norb, j, :);
                    HnumL_tmp(2*i-1, 2*j, :) = H_hr.HnumL(i, j + Norb, :);
                    HnumL_tmp(2*i, 2*j, :) = H_hr.HnumL(i + Norb, j + Norb, :);
                end
            end

        case 'new_to_old'
            % Convert from 'new' format to 'old' format
            for i = 1:Norb
                for j = 1:Norb
                    % Mapping elements from new to old format
                    HnumL_tmp(i, j, :) = H_hr.HnumL(2*i-1, 2*j-1, :);
                    HnumL_tmp(i + Norb, j + Norb, :) = H_hr.HnumL(2*i, 2*j, :);
                end
            end
    end
    
    % Update the Hamiltonian object with the new matrix
    H_hr.HnumL = HnumL_tmp;

    % Optionally write the converted matrix to a file
    if write
        H_hr.Gen_hr('wannier90_hr.converted.dat');
    end
end
