function [mag_ion_seq, sites, Atom_name, Atom_num] = Magion_detect(mag_element_seq_list, Rm, sites, Atom_name, Atom_num, elements)
    % Magion_detect detects and handles magnetic ions (e.g., rare-earth elements) in crystal structure
    % 
    % Input:
    %   mag_element_seq_list - List of atomic numbers for magnetic elements (default: 56:72, Lanthanides)
    %   Rm, sites, Atom_name, Atom_num, elements - Crystal structure parameters (read from POSCAR)
    %
    % Output:
    %   mag_ion_seq - Index of magnetic elements in Atom_name
    %   sites - Reorganized atomic positions (magnetic ions first)
    %   Atom_name - Reorganized atomic names (magnetic ions first)
    %   Atom_num - Reorganized atomic counts (magnetic ions first)

    %% Handle default arguments using inputParser
    p = inputParser;
    addOptional(p, 'mag_element_seq_list', 56:72); % default to Lanthanides
    parse(p, mag_element_seq_list);
    mag_element_seq_list = p.Results.mag_element_seq_list;

    if nargin < 2
        [Rm, sites, Atom_name, Atom_num, elements] = POSCAR_readin('POSCAR', 'vasp');
    end

    %% Create an element query table
    elements.Properties.RowNames = elements.atom_symbol;

    % Display structure information
    fprintf('Crystal contains %d atoms:\n', sum(Atom_num));
    mag_ion_seq = -1;

    %% Loop through all elements and search for magnetic ions
    for i = 1:length(Atom_name)
        tmpdata = table2array(elements(Atom_name{i}, {'atom_number', 'n'}));
        element_seq = tmpdata(1);
        fprintf('    %d %s (Atomic number %d)\n', Atom_num(i), Atom_name{i}, element_seq);
        
        % Check if the element is in the magnetic sequence
        if ismember(element_seq, mag_element_seq_list)
            if mag_ion_seq > 0
                error('Multiple magnetic elements are not supported simultaneously.');
            end
            if Atom_num(i) > 1
                warning('Current version may have limitations for multiple magnetic ions in a unit cell.');
            end
            mag_ion_seq = i; % Record the index of the magnetic element
        end
    end

    %% Validate that a magnetic element was found
    if mag_ion_seq == -1
        error('No magnetic elements detected.');
    else
        fprintf('Magnetic element detected: %s\n', Atom_name{mag_ion_seq});
    end

    %% Generate indices for magnetic atoms
    start_idx = sum(Atom_num(1:mag_ion_seq-1)) + 1;
    end_idx = sum(Atom_num(1:mag_ion_seq));
    mag_indices = start_idx:end_idx;

    %% Separate magnetic component
    Atom_name_mag = Atom_name(mag_ion_seq);
    Atom_num_mag = Atom_num(mag_ion_seq);
    sites_mag = sites(mag_indices);

    % Generate POSCAR for magnetic component
    POSCAR_gen(Rm, sites_mag, Atom_name_mag, Atom_num_mag, 'POSCAR_mag');

    %% Handle non-magnetic component
    Atom_name_others = Atom_name;
    Atom_name_others(mag_ion_seq) = [];
    Atom_num_others = Atom_num;
    Atom_num_others(mag_ion_seq) = [];
    sites_others = sites;
    sites_others(mag_indices) = [];

    % Generate POSCAR for non-magnetic component
    POSCAR_gen(Rm, sites_others, Atom_name_others, Atom_num_others, 'POSCAR_others');

    %% Reorganize structure parameters (magnetic first)
    sites = [sites_mag, sites_others];
    Atom_name = [Atom_name_mag, Atom_name_others];
    Atom_num = [Atom_num_mag, Atom_num_others];

end
