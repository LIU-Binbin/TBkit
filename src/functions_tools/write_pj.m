function num_wan = write_pj(Projector_list)
    % Map for orbital types (n, s, p, d, f, sp, sd, pd, sf, spd, spdf)
    % Corresponding values for orbital types: -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    % Usage: Projector_list will map orbital types to sites

    [~, sites, Atom_name, ~] = POSCAR_read('digits', 16); % Read the POSCAR file

    sites_num = length(sites); % Number of sites in the system
    
    % Map message and orbital type mapping
    maprule_message = ["n s p d f sp sd pd sf spd spdf"; "-1 0 1 2 3 4 5 6 7 8 9"];
    maprule = containers.Map({'n', 's', 'p', 'd', 'f', 'sp', 'sd', 'pd', 'sf', 'spd', 'spdf'}, ...
                             {-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
    
    %% Prompt for Projector List if not provided
    if nargin < 1
        % Display information about each site
        for i = 1:sites_num
            fprintf("%d %4s %f %f %f\n", i, sites(i).name, sites(i).rc1, sites(i).rc2, sites(i).rc3);
        end
        
        Projector_list = -1 * ones(1, sites_num); % Initialize Projector_list to -1 for all sites
        fprintf("%s\n%s\n", maprule_message(1), maprule_message(2)); % Display map rule

        % Prompt user for projectors for each site
        for i = 1:sites_num
            fprintf("%d %4s %f %f %f\n", i, sites(i).name, sites(i).rc1, sites(i).rc2, sites(i).rc3);
            tempchar = input('Please input projector (n s p d f sp sd pd sf spd spdf) of this ion: ', 's');
            Projector_list(i) = maprule(char(tempchar)); % Map input to corresponding orbital type
        end
    end
    
    %% Generate wannier90 projections card (wannier90projector_card)
    maprule2 = containers.Map({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, ...
        {"l=0", "l=1", "l=2", "l=3", "l=0;l=1", "l=0;l=2", "l=1;l=2", "l=0;l=3", "l=0;l=1;l=2", "l=0;l=1;l=2;l=3"});

    fid = fopen('wannier90projector_card', 'w');
    fprintf(fid, "begin projections\n");
    
    % Write projector information for each site
    for i = 1:sites_num
        if Projector_list(i) ~= -1
            fprintf(fid, "f= %f, %f, %f: %s\n", sites(i).rc1, sites(i).rc2, sites(i).rc3, maprule2(Projector_list(i)));
        end
    end
    
    fprintf(fid, "end projections\n");
    fclose(fid);
    
    %% Generate wt_in projectors card (wt_in_projector_card)
    maprule3 = containers.Map({-1, 0, 1, 2, 4, 5, 6, 8}, ...
        {"", "s", "pz px py", "dz2 dxz dyz dx2-y2 dxy", "s pz px py", "s dz2 dxz dyz dx2-y2 dxy", ...
         "pz px py dz2 dxz dyz dx2-y2 dxy", "s pz px py dz2 dxz dyz dx2-y2 dxy"});
    
    maprule4 = containers.Map({-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, ...
        {0, 1, 3, 5, 7, 4, 6, 8, 8, 9, 16});
    
    fid = fopen('wt_in_projector_card', 'w');
    fprintf(fid, "PROJECTORS\n");
    
    num_wan = 0;
    % Write number of Wannier functions for each site
    for i = 1:sites_num
        tempnum = maprule4(Projector_list(i));
        num_wan = num_wan + tempnum; % Update total number of Wannier functions
        fprintf(fid, "%d ", tempnum); 
    end
    fprintf(fid, "\n");

    % Write the projectors for each site
    for i = 1:sites_num
        if Projector_list(i) ~= -1
            fprintf(fid, "%s %s\n", Atom_name(sites(i).nameseq), maprule3(Projector_list(i)));
        end
    end
    
    fclose(fid);
    
    %% Generate wannier90 orbital data (wannier90_orbital.dat)
    digits(16); % Set precision for orbital data output
    maprule5 = containers.Map({-1, 0, 1, 2, 4, 5, 6, 8}, ...
        {[""], ["s"], ["p"], ["d"], ["s", "p"], ["s", "d"], ["p", "d"], ["s", "p", "d"]});
    
    fid = fopen('wannier90_orbital.dat', 'w');
    fprintf(fid, "0 # spinless is 0; spinful is 1\n");
    fprintf(fid, "%d\n", num_wan); % Write the total number of Wannier functions
    
    % Write orbital data for each site
    for i = 1:sites_num
        if Projector_list(i) ~= -1
            wan_stringL = maprule5(Projector_list(i));
            for j = 1:length(wan_stringL)
                fprintf(fid, "%s,", wan_stringL{j});
                fprintf(fid, "%19.16f, %19.16f, %19.16f\n", sites(i).rc1, sites(i).rc2, sites(i).rc3);
            end
        end
    end
    
    fclose(fid);
end
