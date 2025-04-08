function POSCAR_play(sites)
    % POSCAR_play - Prints the site information of a crystal structure
    %
    % Inputs:
    %   sites - An array of structures containing the atomic site information.
    %           Each structure must contain the fields:
    %               'seq' - Sequence index (not used here)
    %               'inseq' - In-sequence index (not used here)
    %               'rc1', 'rc2', 'rc3' - Fractional coordinates of the atom
    %               'name' - The atomic name (element)
    %               'nameseq' - Sequence name of the atom (not used here)
    %
    % Output:
    %   This function will print the following information for each site:
    %   - The site index.
    %   - The atomic element name.
    %   - The fractional coordinates of the atom in the unit cell.

    sites_num = length(sites);  % Number of sites
    fprintf('\nSite Information:\n');
    fprintf('Index  Element   rc1        rc2        rc3\n');
    
    % Loop through all the sites and print their information
    for i = 1:sites_num
        fprintf('%-6d %-8s %-10.6f %-10.6f %-10.6f\n', ...
            i, sites(i).name, sites(i).rc1, sites(i).rc2, sites(i).rc3);
    end
    
    fprintf('\n');  % Add a newline at the end for readability
end
