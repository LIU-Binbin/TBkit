function [Atom_store_smart, nn_store_smart, Rnn, Rnn_map] = nn_smart(Rm, sites, search_range, Accuracy, Rlength_cut)
    % nn_smart calculates the nearest neighbors for a primitive cell.
    % 
    % Inputs:
    %   Rm            - Lattice vectors of the unit cell.
    %   sites         - Atomic sites in the unit cell.
    %   search_range  - Search range for neighbors [x, y, z].
    %   Accuracy      - Precision for distance calculations.
    %   Rlength_cut   - Maximum cutoff distance for neighbors.
    %
    % Outputs:
    %   Atom_store_smart - Struct array containing atomic pair information.
    %   nn_store_smart   - Struct array containing nearest neighbor information.
    %   Rnn              - Unique distances of nearest neighbors.
    %   Rnn_map          - Map of distances to their indices.

    % Set default values for optional inputs if not provided
    arguments
        Rm 
        sites 
        search_range  = [0 0 0];
        Accuracy  = 4;
        Rlength_cut = 15;
    end

    % Initialize variables
    sites_num = size(sites, 2);
    search_rangex = search_range(1);
    search_rangey = search_range(2);
    search_rangez = search_range(3);
    Atom_smart_t = struct('R_fractional_from', [], 'R_fractional_to', [], 'R_fractional_diff', [], 'seq_from', [], 'seq_to', [], 'handyname', []);
    nn_smart_t = struct('seq_from', [], 'seq_to', [], 'nn', []);
    Atom_store_smart = repmat(Atom_smart_t, [sites_num, sites_num]);
    nn_store_smart = repmat(nn_smart_t, [sites_num, sites_num]);
    Rnn_list = [];

    % Generate atomic pair information and nearest neighbor data
    for j = 1:sites_num
        site2 = sites(j); % Home cell
        for i = 1:sites_num
            site1 = sites(i);
            Atom_store_smart(i, j) = Atom_smart_t_gen(site1, site2);
            [Rnn_list_temp, nn_store_smart(i, j)] = nn_smart_t_gen(Atom_store_smart(i, j), Rm, search_rangex, search_rangey, search_rangez, Accuracy, Rlength_cut);
            Rnn_list = [Rnn_list; Rnn_list_temp];
        end
    end

    % Calculate unique distances of nearest neighbors
    Rnn = sort(unique(Rnn_list, 'row'));

    % Create a map of distances to their indices
    Rnn_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for i = 1:length(Rnn)
        Rnn_map(Rnn(i)) = i;
    end

    % Assign levels to nearest neighbors based on their distances
    for j = 1:sites_num
        for i = 1:sites_num
            for k = 1:length(nn_store_smart(i, j).nn)
                nn_store_smart(i, j).nn(k).nn_level = Rnn_map(nn_store_smart(i, j).nn(k).Rlength);
            end
        end
    end
end



function Atom_smart_t = Atom_smart_t_gen(site1, site2)
    % Atom_smart_t_gen generates atomic pair information.
    %
    % Inputs:
    %   site1 - First atomic site.
    %   site2 - Second atomic site.
    %
    % Outputs:
    %   Atom_smart_t - Struct containing atomic pair information.

    Rc1 = [site1.rc1, site1.rc2, site1.rc3];
    Rc2 = [site2.rc1, site2.rc2, site2.rc3];
    Atom_smart_t.R_fractional_from = Rc1;
    Atom_smart_t.R_fractional_to = Rc2;
    Atom_smart_t.R_fractional_diff = -(Rc1 - Rc2);
    Atom_smart_t.seq_from = site1.seq;
    Atom_smart_t.seq_to = site2.seq;
    Atom_smart_t.handyname = strcat(site1.name, ' -> ', site2.name);
end

function [Rnn_list, nn_smart_t] = nn_smart_t_gen(Atom_smart_t, Rm, search_rangex, search_rangey, search_rangez, Accuracy, Rlength_cut)
    % nn_smart_t_gen generates nearest neighbor data for atomic pairs.
    % 
    % Inputs:
    %   Atom_smart_t   - Struct containing atomic pair information (source and destination sites).
    %   Rm             - Lattice vectors of the unit cell.
    %   search_rangex  - Search range in x-direction.
    %   search_rangey  - Search range in y-direction.
    %   search_rangez  - Search range in z-direction.
    %   Accuracy       - Precision for distance calculations.
    %   Rlength_cut    - Maximum cutoff distance for neighbors.
    %
    % Outputs:
    %   Rnn_list       - List of distances of nearest neighbors.
    %   nn_smart_t     - Struct containing nearest neighbor information.
    
    % Initialize the output struct for nearest neighbors
    nn_smart_t.seq_from = Atom_smart_t.seq_from;
    nn_smart_t.seq_to = Atom_smart_t.seq_to;
    
    % Set initial count and calculate the number of possible neighbors
    count = 1;
    reducible_num = (2 * search_rangex + 1) * (2 * search_rangey + 1) * (2 * search_rangez + 1);
    
    % Pre-allocate arrays for storing neighbor data
    Rnn_list = zeros(reducible_num, 1);
    nn_t = struct('R_vector', [], 'R_fractional_diff', [], 'Rlength', [], 'nn_level', []);
    nn = repmat(nn_t, [reducible_num, 1]);
    
    % Loop through possible displacement vectors within the search range
    for Rf_a1 = -search_rangex:search_rangex
        for Rf_a2 = -search_rangey:search_rangey
            for Rf_a3 = -search_rangez:search_rangez
                % Generate the displacement vector
                R_vector = [Rf_a1, Rf_a2, Rf_a3];
                
                % Calculate the length of the displacement vector
                Rlength = norm((R_vector + Atom_smart_t.R_fractional_diff) * Rm);
                
                % Round the calculated distance to the specified accuracy
                Rlength = roundn(Rlength, -Accuracy);
                
                % Check if the distance is within the cutoff range
                if Rlength > 0 && Rlength < Rlength_cut
                    % Store the valid nearest neighbor data
                    nn(count, 1).R_vector = R_vector;
                    nn(count, 1).R_fractional_diff = Atom_smart_t.R_fractional_diff;
                    nn(count, 1).Rlength = Rlength;
                    Rnn_list(count, :) = Rlength;
                    count = count + 1;
                end
            end
        end
    end
    
    % Trim any unused space in the arrays
    if count <= reducible_num
        Rnn_list(count:reducible_num, :) = [];
        nn(count:reducible_num, :) = [];
    end
    
    % Output the nearest neighbor information
    nn_smart_t.nn = nn;
end
