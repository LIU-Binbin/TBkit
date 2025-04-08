function WEIGHTCAR = PROCAR_select(PROCAR_collection, ionslist, orbitalslist, bandsnum, kpointsnum)
    % PROCAR_SELECT Computes the weighted sum of projected densities of states (pDOS).
    %
    %   WEIGHTCAR = PROCAR_SELECT(PROCAR_collection, ionslist, orbitalslist, bandsnum, kpointsnum)
    %
    % Inputs:
    %   PROCAR_collection - Struct array containing PROCAR data for each ion and orbital.
    %   ionslist          - Vector of ion indices to include in the calculation.
    %   orbitalslist      - Vector of orbital indices to include in the calculation.
    %   bandsnum          - Total number of bands.
    %   kpointsnum        - Total number of k-points.
    %
    % Outputs:
    %   WEIGHTCAR         - Matrix of weighted pDOS with dimensions [bandsnum, kpointsnum].
    %
    % Note:
    %   Ensure that PROCAR_collection is a struct array with fields corresponding to ions and orbitals.
    %   ionslist and orbitalslist should be vectors of indices referring to the desired ions and orbitals.
    %   bandsnum and kpointsnum should be scalar integers representing the total number of bands and k-points.
    %
    % Author: [Your Name]
    % Date: [Date]
    % Version: 1.0
    % License: [License Information]
    % Contact: [Your Contact Information]

    % Validate input parameters
    if nargin < 5
        error('Insufficient input arguments. Please provide PROCAR_collection, ionslist, orbitalslist, bandsnum, and kpointsnum.');
    end

    % Initialize the WEIGHTCAR matrix with zeros
    WEIGHTCAR = zeros(bandsnum, kpointsnum);

    % Loop over each ion in ionslist
    for i = 1:length(ionslist)
        % Loop over each orbital in orbitalslist
        for j = 1:length(orbitalslist)
            % Add the weighted pDOS for the current ion and orbital to the total
            WEIGHTCAR = WEIGHTCAR + PROCAR_collection(ionslist(i), orbitalslist(j)).WEIGHTCAR;
        end
    end
end
