function orbital_out = nanowire_orb(H_hr,fin_dir,vacuum_mode,options)
% NANOWIRE_ORB Generate orbital positions for nanowire structure
%
%   ORBITAL_OUT = NANOWIRE_ORB(H_HR,FIN_DIR,VACUUM_MODE,OPTIONS) generates
%   orbital positions for nanowire construction
%
%   Inputs:
%       H_hr - HR object containing initial orbital information
%       fin_dir - Fin direction vector [default: [10 10 1]]
%       vacuum_mode - Boolean flag for vacuum mode [default: false]
%       options.fast - Boolean flag for fast mode [default: true]
%   Output:
%       orbital_out - Generated orbital positions
%
%   Notes:
%       - Handles both vacuum and non-vacuum modes
%       - Supports supercell generation in non-fast mode
%       - Processes orbital positions along specified fin directions
arguments
H_hr HR;
fin_dir = [10 10 1];
vacuum_mode = false;
options.fast = true;
end
if isempty(H_hr.orbL)
orbital_init = zeros(H_hr.WAN_NUM,3);
else
orbital_init = H_hr.orbL;
end
if ~vacuum_mode
orbital_out = orbital_init;
for i = 1:H_hr.Dim
Nslab = fin_dir(i);
if  Nslab == 0
Nslab = 1;
end
count = 0;
WAN_NUM_tmp = size(orbital_out,1);
fin_orb = zeros(WAN_NUM_tmp*Nslab,3);
for inum = 1:Nslab
for j = 1:WAN_NUM_tmp
count =count +1;
orb_tmp=orbital_out(j,:);
orb_tmp(i)= (orb_tmp(i) + double(inum-1)) / Nslab;
fin_orb(count,:)=orb_tmp;
end
end
orbital_out = fin_orb;
end
Ns = [1 0 0;0 1 0;0 0 1];
Ns = Ns.*fin_dir;
fin_dir_list = [0 0 0];
if ~options.fast
[Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp]=HR.POSCAR_read();
H_hr.supercell(Ns,'POSCAR_super_fin',Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp,fin_dir_list);
end
else
orbital_out = orbital_init;
for i = 1:H_hr.Dim
Nslab = fin_dir(i);
if  Nslab == 0
Nslab = 1;
end
count = 0;
WAN_NUM_tmp = size(orbital_out,1);
fin_orb = zeros(WAN_NUM_tmp*Nslab,3);
for inum = 1:Nslab
for j = 1:WAN_NUM_tmp
count =count +1;
orb_tmp=orbital_out(j,:);
orb_tmp(i)= (orb_tmp(i) + double(inum-1))/Nslab;
fin_orb(count,:)=orb_tmp;
end
end
orbital_out = fin_orb;
end
Ns = [1 0 0;0 1 0;0 0 1];
Ns = Ns.*fin_dir;
fin_dir_list = double(fin_dir>1);
if ~options.fast
[Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp]=HR.POSCAR_read();
H_hr.supercell(Ns,'POSCAR_super_fin',Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp,fin_dir_list);
else
[Rm_tmp,sites_tmp,Atom_name_tmp,Atom_num_tmp]=POSCAR_read();
end
Rm_tmp = Ns*Rm_tmp;
Rmlength1 = norm (Rm_tmp(1,:));
Rmlength2 = norm (Rm_tmp(2,:));
Rmlength3 = norm (Rm_tmp(3,:));
Rm_s_fin_add = [10*Rm_tmp(1,:)*fin_dir_list(1)/Rmlength1;...
10*Rm_tmp(2,:)*fin_dir_list(2)/Rmlength2;...
10*Rm_tmp(3,:)*fin_dir_list(3)/Rmlength3];
Rm_s_fin = Rm_tmp + Rm_s_fin_add ;
Rc_s_fin_add = [1/2, 1/2 ,1/2] ;
Rr_s_fin_add = Rc_s_fin_add * Rm_s_fin_add;
[nfinorb,~ ]= size(fin_orb);
for  i = 1:nfinorb
Rr_orb = fin_orb(i,:)*Rm_tmp;
Rr_s_fin = Rr_orb + Rr_s_fin_add;
Rc_s_fin = Rr_s_fin / Rm_s_fin;
fin_orb(i,:) = Rc_s_fin ;
end
orbital_out = fin_orb;
end
end
