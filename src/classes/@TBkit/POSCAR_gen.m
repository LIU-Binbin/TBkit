function [Rm,sites,Atom_name,Atom_num] = POSCAR_gen(TBkitobj,filename,Rm,sites,Atom_name,Atom_num)
%POSCAR_GEN Generate VASP POSCAR file from TBkit object
%   [RM, SITES, ATOM_NAME, ATOM_NUM] = POSCAR_GEN(TBKITOBJ, FILENAME, ...)
%   generates a VASP-format POSCAR file from a TBkit object.
%
%   Inputs:
%       TBkitobj   - TBkit object containing structural information
%       filename   - Output filename (default: 'POSCAR')
%       Rm         - Lattice vectors (default: from TBkitobj)
%       sites      - Atomic sites (default: from TBkitobj)
%       Atom_name  - Atom names (default: from TBkitobj)
%       Atom_num   - Atom numbers (default: from TBkitobj)
%
%   Outputs:
%       Rm         - Lattice vectors
%       sites      - Atomic sites
%       Atom_name  - Atom names
%       Atom_num   - Atom counts
%
%   Note:
%       Writes output in VASP POSCAR format with direct coordinates
arguments
    TBkitobj
    filename = 'POSCAR';
    Rm = TBkitobj.Rm;
    sites = TBkitobj.sites;
    Atom_name = TBkitobj.Atom_name;
    Atom_num = TBkitobj.Atom_num;
end
title = "POSCAR Gen by vasplib";
%% write POSCAR
% Initialize variables
%filename = "POSCAR_"+num2str(term2)+".vasp";
fileID = fopen(filename,'w');
%% Rm
a_crystal_constance=1;
%%
fprintf(fileID,"%s\n",title);
fprintf(fileID,"%d\n",a_crystal_constance);
%fprintf(fileID,"  ",Rm(i,j));
for i=1:3
    for j=1:3
        fprintf(fileID,"  %f",Rm(i,j));
    end
    fprintf(fileID,"\n");
end
for i=1:length(Atom_name)
    fprintf(fileID,"%s ",Atom_name(i));
end
fprintf(fileID,"\n  ");
for i=1:length(Atom_num)
    fprintf(fileID,"%d ",Atom_num(i));
end
fprintf(fileID,"\n");
fprintf(fileID,"Direct\n  ");
% sites
[~,sites_num]=size(sites);
for i=1:sites_num
    fprintf(fileID,"%f  ",mod(sites(i).rc1,1));
    fprintf(fileID,"%f  ",mod(sites(i).rc2,1));
    fprintf(fileID,"%f  ",mod(sites(i).rc3,1));
    fprintf(fileID,"\n  ");
end
fclose(fileID);
end

