function test = PBAND_DAT_gen(outputfile, EIGENCAR, WEIGHTCAR, klist_l)
    % PBAND_DAT_GEN Generates a projected band structure data file.
    %
    %   test = PBAND_DAT_GEN(outputfile, EIGENCAR, WEIGHTCAR, klist_l)
    %
    % Inputs:
    %   outputfile - Name of the output data file.
    %   EIGENCAR   - Eigenvalue matrix from VASP calculations.
    %   WEIGHTCAR  - Corresponding weight matrix.
    %   klist_l    - List of k-points.
    %
    % Outputs:
    %   test       - Status code (0 indicates success).
    %
    % Note:
    %   Ensure that EIGENCAR, WEIGHTCAR, and klist_l are loaded before calling this function.
    %
    % Author: [Your Name]
    % Date: [Date]
    % Version: 1.0
    % License: [License Information]
    % Contact: [Your Contact Information]

    % Validate input parameters
    if nargin < 4
        error('Insufficient input arguments. Please provide outputfile, EIGENCAR, WEIGHTCAR, and klist_l.');
    end

    % Get the number of bands and k-points
    [bandsnum, kpointsnum] = size(EIGENCAR);

    % Open the output file for writing
    try
        PBANDDATAFILE = fopen(outputfile, 'w');
        if PBANDDATAFILE == -1
            error('Failed to open file: %s', outputfile);
        end
    catch
        error('File operation failed: %s', outputfile);
    end

    % Write data to the file
    for i = 1:bandsnum
        fprintf(PBANDDATAFILE, '#%9s %10s %10s   ------BAND: %d\n', 'KPOINTS', 'EIGEN', 'WEIGHT', i);
        for j = 1:kpointsnum
            k_point = klist_l(j);  % Avoid redundant computation
            fprintf(PBANDDATAFILE, '%10.6f %10.6f %10.6f\n', k_point, EIGENCAR(i, j), WEIGHTCAR(i, j));
        end
    end

    % Close the file
    fclose(PBANDDATAFILE);

    % Return status code
    test = 0;
end
