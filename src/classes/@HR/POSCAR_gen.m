function [H_hr_bk] = POSCAR_gen(H_hr,filename,options)
% POSCAR_GEN Generate VASP POSCAR file from HR object
%
%   [H_HR_BK] = POSCAR_GEN(H_HR,FILENAME,OPTIONS) generates a VASP POSCAR
%   file from orbital information
%
%   Inputs:
%       H_hr - HR object containing orbital data
%       filename - Output filename [default: 'POSCAR_gen.vasp']
%       options.title - File header [default: includes date]
%       options.a_crystal_constance - Scaling factor [default: 1]
%       options.vacuum - Flag to add vacuum [default: false]
%       options.fin_dir_list - Fin directions [default: [1 1 0]]
%       options.vacuum_length - Vacuum thickness [default: 10]
%   Output:
%       H_hr_bk - Backup of original HR object
%
%   Notes:
%       - Handles vacuum layer addition
%       - Reads element information from 'elements.txt'
%       - Preserves original HR object
arguments
H_hr HR;
filename = 'POSCAR_gen.vasp';
options.title = "PARCHG Gen by TBkit. " + string(date);
options.a_crystal_constance = 1;
options.vacuum = false;
options.fin_dir_list = [1 1 0];
options.vacuum_length = 10;
end
if options.vacuum
[H_hr.orbL,H_hr.Rm] = TBkit.AddVacuumLayer(H_hr.orbL,H_hr.Rm,options.fin_dir_list,'vacuum_length',options.vacuum_length);
end
H_hr_bk = H_hr;
quantumList = H_hr.quantumL;
H_hr = H_hr.reseq(quantumList(:,end)>=0 |isnan(quantumList(:,end)));
elementList = H_hr.elementL;
[elementList,sortL] = sort(elementList);
H_hr = H_hr.reseq(sortL);
orbList = H_hr.orbL;
try
elements = readtable('elements.txt');
catch
elements = [];
end
if ~isempty(elements)
end
unique_elementL = unique(elementList);
AtomTypes = length(unique_elementL);
Atom_num = zeros(1,AtomTypes);
Atom_name = "";
Atom_name = repmat(Atom_name,size(Atom_num));
Rm = H_hr.Rm;
for i = 1:AtomTypes
atom_number = unique_elementL(i);
Atom_num(i) = sum(elementList == atom_number);
rows = elements.atom_number ==atom_number;
Atom_name(i) = string((elements.atom_symbol(rows)));
end
fileID = fopen(filename,'w');
fprintf(fileID,"%s\n",options.title);
fprintf(fileID,"%d\n",options.a_crystal_constance);
%fprintf(fileID,"  ",Rm(i,j));
for i=1:H_hr.Dim
for j=1:H_hr.Dim
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
for i=1:size(orbList,1)
fprintf(fileID,"%f  ",mod(orbList(i,1),1));
fprintf(fileID,"%f  ",mod(orbList(i,2),1));
fprintf(fileID,"%f  ",mod(orbList(i,3),1));
%             %fprintf(fileID,"%s\n  ",sites(i).name);
%             fprintf(fileID,"\n  ");
%             fprintf(fileID,"\n  ");
fprintf(fileID,"\n  ");
end
fclose(fileID);
[H_hr_bk.Rm,H_hr_bk.sites,H_hr_bk.Atom_name,H_hr_bk.Atom_num,~]=POSCAR_read(filename);
end
