function [Rm,sites,Atom_name,Atom_num,elements]=POSCAR_read(filename,mode,options)
%POSCAR_READ Parse crystal structure information from POSCAR files
%
%   Syntax:
%       [Rm,sites,Atom_name,Atom_num,elements] = POSCAR_read(filename,mode,options)
%
%   Description:
%       Reads VASP POSCAR files and extracts crystal structure information
%       including lattice vectors, atomic positions, and orbital data.
%
%   Inputs:
%       filename - Path to POSCAR file (default='POSCAR')
%       mode     - Parsing mode ('vasp','symble','tbsk','tbsym','sym','list')
%       options  - Structure with parsing options (digits)
%
%   Outputs:
%       Rm        - 3×3 lattice vectors matrix
%       sites     - Structure array of atomic positions
%       Atom_name - Cell array of atomic symbols
%       Atom_num  - Array of atom counts
%       elements  - Periodic table data
%
%   See also: POSCAR_cell_read, POSCAR_gen
%
% get poscar information from vasp and others
% * Label:
%
% Description of the Function:
%
arguments
    filename string ='POSCAR';
    mode char {mustBeMember(mode,{'symble','vasp','tbsk','tbsym','sym','list'})} = 'vasp';
    options.digits = 15;
end
%--------  init  --------
try
    elements = readtable('elements.txt');
catch
    elements = [];
end
%--------  chek  --------
if ~exist(filename,'file')
    fprintf('No such file: %s\n',filename);
    error('file not exist!');
end
%--------  juge  --------
switch mode
    case {'vasp','symble'}
        site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);
        formatSpec = '%s%s%s%s%s%[^\n\r]';
    case {'tbsk','tbsym'}
        site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[],'element',[],'orb',[],'orb_sym',[],'spin',[]);
        formatSpec = '%s%s%s%s%s%s%s%[^\n\r]';
    case 'sym'
        site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[],'element',[],'orb',[],'orb_sym',[],'spin',[],'Jz',[]);
        formatSpec = '%s%s%s%s%s%s%s%s%[^\n\r]';
    case 'list'
        site=struct('seq',[],'rc1',[],'rc2',[],'rc3',[],'Hue',[],'surf_level',[],'hing_level',[]);

end
%--------  init  --------
digits(options.digits); %take care
POSCAR = TBkit.POSCAR_cell_read(filename,formatSpec);
% a_crystal_constance=str2double(char(POSCAR(1,1)));
% Rm
a1=[str2double(char(POSCAR(2,1))) str2double(char(POSCAR(2,2))) str2double(char(POSCAR(2,3)))];
a2=[str2double(char(POSCAR(3,1))) str2double(char(POSCAR(3,2))) str2double(char(POSCAR(3,3)))];
a3=[str2double(char(POSCAR(4,1))) str2double(char(POSCAR(4,2))) str2double(char(POSCAR(4,3)))];
Rm=[a1;a2;a3];
%
%Rm = (Rm)/det(Rm);
% read atom
n_atom=length(POSCAR(5,:));
for i=1:n_atom
    if POSCAR(5,i) ~= ""
        Atom_name(i) = POSCAR(5,i);
    end
end
% need
for i=1:length(Atom_name)
    Atom_num(i)=str2double(char(POSCAR(6,i)));
end
%site_num
sites_num=sum(Atom_num);
sites=repmat(site,[1 sites_num]);
% coordinate pattern
% Coordinates_pattern=POSCAR(7,1);
% temp : frac support
% struct
sequence=1;
n_atom=length(Atom_name);
labelcut_list = TBkit.labelcut_list_gen(Atom_num);

%first
if n_atom >= 1
    for i=1:n_atom
        inseq =1;
        for j=labelcut_list(i,1):labelcut_list(i,2)
            if j > size(POSCAR,1)
                error('Your POSCAR 7th line mismatches the fractional coordinates of each orbital  ')
            end
            %id
            sites(sequence).seq=sequence;
            %inid
            sites(sequence).inseq=inseq ;
            %name
            sites(sequence).nameseq=i;
            sites(sequence).name=Atom_name(i)+num2str(sites(sequence).inseq);
            sites(sequence).ion_num=Atom_num(i);
            %sites(sequence).ion_type_num=i;
            %incoordinate
            sites(sequence).rc1=str2double(char(POSCAR(j,1)));
            sites(sequence).rc2=str2double(char(POSCAR(j,2)));
            sites(sequence).rc3=str2double(char(POSCAR(j,3)));
            if strcmp(mode,'tbsk')
                sites(sequence).element = string(POSCAR(j,4));
                sites(sequence).orb = string(POSCAR(j,5));
                sites(sequence).orb_sym = string(POSCAR(j,6));
                if length(POSCAR(j,:)) > 6
                    sites(sequence).spin = str2double(POSCAR(j,7));
                else
                    sites(sequence).spin = 1;
                end
            elseif strcmp(mode,'tbsym')
                sites(sequence).element = string(POSCAR(j,4));
                sites(sequence).orb = string(POSCAR(j,5));
                sites(sequence).orb_sym = string(POSCAR(j,6));
                sites(sequence).spin = str2double(POSCAR(j,7));
            elseif strcmp(mode,'sym')
                sites(sequence).element = string(POSCAR(j,4));
                sites(sequence).orb = string(POSCAR(j,5));
                sites(sequence).orb_sym = string(POSCAR(j,6));
                sites(sequence).spin = str2double(POSCAR(j,7));
                sites(sequence).J = str2double(POSCAR(j,8));
            end
            %
            sequence=sequence+1;
            inseq =inseq + 1;

        end
    end
end

if length(POSCAR)>j && strcmp(mode,'mag')
    %disp('mag_mode');
    sequence=1;
    beginline=j+1;
    %first
    for j=beginline:beginline+Atom_num(1)-1
        %id
        sites(sequence).mag=double(POSCAR(j,1:3));
        sequence=sequence+1;
    end
    %other
    if n_atom >= 2
        for i=2:n_atom
            beginline=beginline+Atom_num(i-1);
            for j=beginline:beginline+Atom_num(i)-1
                %id
                sites(sequence).mag=double(POSCAR(j,1:3));
                sequence=sequence+1;
            end
        end
    end
end
end