function [Rm, sites, Atom_name, Atom_num, elements, a_crystal_constance] = POSCAR_read(filename, mode, options)
%POSCAR_READ Reads crystal structure information from POSCAR-type files
%   This function parses VASP-compatible structure files and related formats,
%   extracting lattice vectors, atomic positions, and additional orbital/spin
%   information depending on the specified format mode.
%
%   Inputs:
%       filename  - Path to input file (default: 'POSCAR')
%       mode      - File format mode: 'vasp', 'tbsk', 'tbsym', 'list' (default: 'vasp')
%       options   - Struct with processing options:
%                   digits: Precision for coordinate values (default: 6)
%
%   Outputs:
%       Rm                  - 3x3 matrix of lattice vectors [a1; a2; a3]
%       sites               - Structure array containing atomic positions and metadata
%       Atom_name           - Cell array of element names for each atom type
%       Atom_num            - Array containing number of atoms per element
%       elements            - Periodic table data from elements.txt
%       a_crystal_constance - Scaling factor for lattice vectors

%% Argument validation and initialization
arguments
    filename = 'POSCAR';
    mode = 'vasp';
    options.digits = 6;
end

% Load periodic table data for element lookup
elements = readtable('elements.txt');
digits(options.digits); % Set numeric precision for coordinates

%% File existence check
if ~exist(filename,'file')
    fprintf('No such file: %s\n',filename);
    error('POSCAR_READ:FileNotFound', 'Target file does not exist');
end

%% Initialize data structures based on file format
switch mode
    case 'vasp'
        site = struct('seq',[], 'inseq',[], 'rc1',[], 'rc2',[], 'rc3',[],...
            'name',[], 'nameseq',[]);
        formatSpec = '%s%s%s%s%s%[^\n\r]';
        
    case 'tbsk'
        site = struct('seq',[], 'inseq',[], 'rc1',[], 'rc2',[], 'rc3',[],...
            'name',[], 'nameseq',[], 'orb',[], 'orb_sym',[]);
        formatSpec = '%s%s%s%s%s%s%[^\n\r]';
        
    case 'tbsym'
        site = struct('seq',[], 'inseq',[], 'rc1',[], 'rc2',[], 'rc3',[],...
            'name',[], 'nameseq',[], 'element',[], 'orb',[], 'orb_sym',[], 'spin',[]);
        formatSpec = '%s%s%s%s%s%s%s%[^\n\r]';
        
    case 'list'
        site = struct('seq',[], 'rc1',[], 'rc2',[], 'rc3',[],...
            'Hue',[], 'surf_level',[], 'hing_level',[]);
end

%% File parsing core
POSCAR = POSCAR_cell_read(filename, formatSpec);

%% Lattice parameters extraction
a_crystal_constance = str2double(char(POSCAR(1,1))); % Scaling factor
Rm = [str2double((POSCAR(2,1:3)));    % Lattice vector a1
       str2double((POSCAR(3,1:3)));    % Lattice vector a2
       str2double((POSCAR(4,1:3)))];   % Lattice vector a3

%% Atomic species processing
Atom_name = POSCAR(5, ~cellfun(@isempty, POSCAR(5,:))); % Extract element names
Atom_num = cellfun(@str2double, POSCAR(6,1:length(Atom_name))); % Atom counts

%% Atomic position data organization
sites_num = sum(Atom_num); % Total atomic sites
sites = repmat(site, [1, sites_num]); % Preallocate structure array
labelcut_list = labelcut_list_gen(Atom_num); % Generate line indices

%% Populate atomic site information
sequence = 1; % Global site index counter
for i = 1:length(Atom_name)
    inseq = 1; % Per-element site index
    for j = labelcut_list(i,1):labelcut_list(i,2)
        sites(sequence).seq = sequence;
        sites(sequence).inseq = inseq;
        sites(sequence).nameseq = i;
        sites(sequence).name = sprintf('%s%d', Atom_name(i), inseq);
        sites(sequence).rc1 = str2double((POSCAR(j,1)));
        sites(sequence).rc2 = str2double((POSCAR(j,2)));
        sites(sequence).rc3 = str2double((POSCAR(j,3)));
        
        % Handle extended format data
        if strcmp(mode, 'tbsk')
            sites(sequence).orb = string(POSCAR(j,5));
            sites(sequence).orb_sym = sym(POSCAR(j,6));
        elseif strcmp(mode, 'tbsym')
            sites(sequence).element = string(POSCAR(j,4));
            sites(sequence).orb = string(POSCAR(j,5));
            sites(sequence).orb_sym = string(POSCAR(j,6));
            sites(sequence).spin = str2double(POSCAR(j,7));
        end
        
        sequence = sequence + 1;
        inseq = inseq + 1;
    end
end

%% Magnetic moment handling (if present)
if length(POSCAR) > j && strcmp(mode, 'mag')
    sequence = 1;
    beginline = j + 1;
    % ... (magnetic moment parsing logic here)
end
end

%% Helper function: Generate line index ranges for atomic positions
function labelcut_list = labelcut_list_gen(Atom_num)
%LABELCUT_LIST_GEN Generates line index ranges for atomic position blocks
    n = length(Atom_num);
    beginline = 8; % POSCAR format standard starting line
    sum_n = beginline;
    
    for i = 1:n
        sum_n2 = sum_n + Atom_num(i) - 1;
        labelcut_list(i,:) = [sum_n, sum_n2];
        sum_n = sum_n2 + 1;
    end
end

%% Helper function: Core file reader
function POSCAR = POSCAR_cell_read(filename, formatSpec)
%POSCAR_CELL_READ Low-level text parser for POSCAR-formatted files
    startRow = 2;
    delimiter = {'\t',' '};
    
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,...
        'MultipleDelimsAsOne', true, 'TextType', 'string',...
        'HeaderLines', startRow-1, 'ReturnOnError', false,...
        'EndOfLine', '\r\n');
    fclose(fileID);
    
    % Convert to uniform string array
    raw = repmat({''}, length(dataArray{1}), length(dataArray)-1);
    for col = 1:length(dataArray)-1
        raw(1:length(dataArray{col}), col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
    end
    POSCAR = string(raw(:,:));
end
