function [Rm, sites] = POSCAR_rmRc(rc_rm_list, Accuracy, filename, Rm, sites, Atom_name, Atom_num)
    % POSCAR_rmRc - Removes specific atomic sites based on fractional coordinates.
    %
    % Inputs:
    %   rc_rm_list   - List of fractional coordinates to remove (N x 3 matrix).
    %   Accuracy     - Precision level for comparing coordinates (default: 4).
    %   filename     - Output filename for the modified POSCAR (default: 'POSCAR_rmRc').
    %   Rm           - Current crystal matrix.
    %   sites        - Structure array containing atomic site information.
    %   Atom_name    - Cell array of atom names.
    %   Atom_num     - Array of atom counts for each element.
    %
    % Outputs:
    %   Rm           - Updated crystal matrix.
    %   sites        - Updated structure array after removal.

    %% Check input arguments
    if nargin < 4
        POSCAR_read;  % Read POSCAR if not provided
    end
    if nargin < 3
        filename = 'POSCAR_rmRc';  % Default filename
    end
    if nargin < 2
        Accuracy = 4;  % Default accuracy level
    end

    % Extract fractional coordinates from sites structure
    Rc_list = [[sites.rc1]', [sites.rc2]', [sites.rc3]'];

    %% Round coordinates to match desired accuracy
    rc_rm_list = round(rc_rm_list .* 10^Accuracy) / 10^Accuracy;
    rc_rm_list = unique(rc_rm_list, 'rows');  % Ensure no duplicate coordinates

    Rc_list = round(Rc_list .* 10^Accuracy) / 10^Accuracy;  % Round site coordinates

    %% Find matching coordinates to remove
    seq_list = find_unique_list_in_a_list(rc_rm_list, Rc_list);
    
    %% Handle atom type list
    atomtype_list = {sites(seq_list).nameseq}';  % Extract atom types for the selected sites

    % Create atom type list
    Atom_type_list = 1:length(Atom_name);  % Atom type index list
    atomtype_list = [Atom_type_list; atomtype_list];

    %% Remove duplicate atoms from Atom_num
    repeat_list = histc(atomtype_list, unique(atomtype_list)) - 1;
    Atom_num = Atom_num - repeat_list';  % Adjust atom counts based on removals

    %% Remove empty atom types (if any)
    [atom_unique, atom_seq_list] = ismember([0], Atom_num', 'rows');
    if sum(atom_unique) > 0
        Atom_num(atom_seq_list) = [];
        Atom_name(atom_seq_list) = [];
    end

    %% Remove selected atomic sites
    sites(seq_list) = [];

    %% Regenerate the POSCAR file with updated atomic site information
    [Rm, sites] = POSCAR_gen(Rm, sites, Atom_name, Atom_num, filename);
end

%% Helper Function: Find unique entries from list_a in list_b
function seq_list = find_unique_list_in_a_list(list_a, list_b)
    % This function finds matching rows of list_a within list_b.
    % Returns the indices of rows in list_b that match rows in list_a.
    
    seq_list = [];
    
    if isempty(list_a)
        disp('No coordinates to remove');
        return;
    end
    
    % Find the indices of matching rows in list_b
    [all_one, seq_list_init] = ismember(list_a, list_b, 'rows');
    
    % Extract only the matching sequences
    nseq_list = length(all_one);
    for i = 1:nseq_list
        if all_one(i) == 1
            seq_list = [seq_list; seq_list_init(i, :)];
        end
    end
end
