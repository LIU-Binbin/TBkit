function [POSCAR_file, POSCAR_syminfor_file, ID] = SymInfor_gen(POSCAR_name, phonopy_run)
    % SymInfor_gen - Generate symmetry information using phonopy.
    %
    % Syntax:
    %   [POSCAR_file, POSCAR_syminfor_file, ID] = SymInfor_gen(POSCAR_name, phonopy_run)
    %
    % Input:
    %   POSCAR_name - The name of the POSCAR file (e.g., 'POSCAR')
    %   phonopy_run - Path to the phonopy executable (optional, default set to a specific path)
    %
    % Output:
    %   POSCAR_file         - The original POSCAR file name
    %   POSCAR_syminfor_file - The generated symmetry information file name
    %   ID                   - The ID extracted from the POSCAR file name
    
    % Default phonopy executable path if not provided
    if nargin < 2
        phonopy_run = '/Users/parkman/Documents/TOOLs/miniconda2/envs/my_pymatgen/bin/phonopy';
    end
    
    % Generate symmetry information using phonopy
    cmd = sprintf('!%s --tolerance 0.01 --symmetry -c %s > SymInfor_%s', phonopy_run, POSCAR_name, POSCAR_name);
    eval(cmd);  % Run phonopy command to get symmetry info
    
    % Try to handle the POSCAR file and move to the current directory
    try
        copyfile('PPOSCAR', 'POSCAR');
    catch
        disp('Error: Invalid POSCAR file');
        ID = -1;
        POSCAR_file = '';
        POSCAR_syminfor_file = '';
        return;
    end
    
    % Run phonopy symmetry analysis on the 'POSCAR' file
    cmd = sprintf('!%s --tolerance 0.01 --symmetry -c POSCAR > SymInfor_%s', phonopy_run, POSCAR_name);
    eval(cmd);  % Run the symmetry analysis
    
    % Clean up temporary files
    delete('BPOSCAR');
    delete('POSCAR');
    
    % Move 'PPOSCAR' back to the original POSCAR file
    try
        movefile('PPOSCAR', POSCAR_name);
    catch
        disp('Error: Unable to restore the original POSCAR file');
        return;
    end
    
    % Set output variables
    POSCAR_file = POSCAR_name;
    POSCAR_syminfor_file = sprintf('SymInfor_%s', POSCAR_name);
    ID = str2double(strrep(POSCAR_file, 'POSCAR.', ''));
end
