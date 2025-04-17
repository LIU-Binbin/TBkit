function Hnum_list_wire_iz =  Hnum_list_wire_iz_gen(Hnum_list_xy_iz,vertor_list_xy_iz,iz,WAN_NUM,WAN_NUM_x,WAN_NUM_y,Nslab,Type)
% HNUM_LIST_WIRE_IZ_GEN Generate wire Hamiltonian block for specific z-index
%
%   Hnum_list_wire_iz = HNUM_LIST_WIRE_IZ_GEN(Hnum_list_xy_iz,vertor_list_xy_iz,iz,WAN_NUM,WAN_NUM_x,WAN_NUM_y,Nslab,Type)
%   constructs a single block of the nanowire Hamiltonian for a given z-index.
%
%   INPUT ARGUMENTS:
%       Hnum_list_xy_iz - Hamiltonian blocks for xy plane
%       vertor_list_xy_iz - Vector list for xy plane
%       iz - Current z-index
%       WAN_NUM - Number of Wannier orbitals
%       WAN_NUM_x - Number of orbitals along x
%       WAN_NUM_y - Number of orbitals along y
%       Nslab - Supercell dimensions
%       Type - Storage type ('sparse' or 'mat')
%
%   OUTPUT ARGUMENTS:
%       Hnum_list_wire_iz - Generated Hamiltonian block
%
%   NOTES:
%       - Handles both sparse and dense matrix formats
%       - Includes progress reporting for dense matrices
%
%   SEE ALSO:
%       HR, Hnanowire_gen
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]
Hnum_list_wire_iz = zeros(WAN_NUM_y,WAN_NUM_y);
[NRPTS_xy,~] = size(vertor_list_xy_iz);
[unique_y,unique_label_y]= unique(vertor_list_xy_iz(:,2),'rows');
cutlist_y = HR.unique_label2cutlist(unique_label_y,NRPTS_xy);
NRPTS_y = length(unique_y);
for iy = 1:NRPTS_y
fprintf('%d th NRPT z ---- Gen (%d/%d) NRPT y \n',iz,iy,NRPTS_y);
vertor_list_x_iy_iz = vertor_list_xy_iz(cutlist_y(iy,1):cutlist_y(iy,2),:);
if strcmp(Type,'sparse')
Hnum_list_x_iy_iz   = Hnum_list_xy_iz(cutlist_y(iy,1):cutlist_y(iy,2));
else
Hnum_list_x_iy_iz   = Hnum_list_xy_iz(:,:,cutlist_y(iy,1):cutlist_y(iy,2));
end
[NRPTS_x,~] = size(vertor_list_x_iy_iz);
Hnum_list_iy_iz =zeros(WAN_NUM_x,WAN_NUM_x);
if strcmp(Type,'sparse')
for ix = 1:NRPTS_x
[ilist,jlist,amplist] = find(Hnum_list_x_iy_iz{ix});
nhopping = length(amplist);
fprintf('%d th NRPT z ---- %d th NRPT y ---- Gen (%d/%d) NRPT x \n',iz,iy,ix,NRPTS_x);
ind_R_x = vertor_list_x_iy_iz(ix,:);
jump_fin_x=ind_R_x(1);
for ih = 1:nhopping
i = ilist(ih);
j = jlist(ih);
amp = amplist(ih);
for icur_sc_vec = 1:Nslab(1)
hi= i + (icur_sc_vec-1)*WAN_NUM ;
hj= j + (icur_sc_vec+jump_fin_x-1)*WAN_NUM ;
to_add=1;
if hj <= 0 || hj > WAN_NUM_x
to_add=0;
end
if to_add == 1
Hnum_list_iy_iz(hi,hj) = amp;
end
end
end
end
else
pb = CmdLineProgressBar('Gen NRPT x: ');
for ix = 1:NRPTS_x
pb.print(ix,NRPTS_x);
ind_R_x = vertor_list_x_iy_iz(ix,:);
jump_fin_x=ind_R_x(1);
for  i = 1:WAN_NUM
for j = 1:WAN_NUM
amp = Hnum_list_x_iy_iz(i,j,ix);
if norm(amp) > 0
for icur_sc_vec = 1:Nslab(1)
hi= i + (icur_sc_vec-1)*WAN_NUM ;
hj= j + (icur_sc_vec+jump_fin_x-1)*WAN_NUM ;
to_add=1;
if hj <= 0 || hj > WAN_NUM_x
to_add=0;
end
if to_add == 1
Hnum_list_iy_iz(hi,hj) = amp;
end
end
end
end
end
end
pb.delete
end
jump_fin_y =  unique_y(iy);
[ilist,jlist,amplist] = find(Hnum_list_iy_iz);
nhopping = length(amplist);
for ih = 1:nhopping
i = ilist(ih);
j = jlist(ih);
amp = amplist(ih);
for icur_sc_vec = 1:Nslab(2)
hi= i + (icur_sc_vec-1)*WAN_NUM_x ;
hj= j + (icur_sc_vec+jump_fin_y-1)*WAN_NUM_x ;
to_add=1;
if hj <= 0 || hj > WAN_NUM_y
to_add=0;
end
if to_add == 1
Hnum_list_wire_iz(hi,hj) =  amp;
end
end
end
end
Hnum_list_wire_iz = sparse(Hnum_list_wire_iz);
end
