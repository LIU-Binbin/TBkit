function H_hr = Hnanowire_gen(H_hr,Nslab,np,vacuum_mode,options)
% HNANOWIRE_GEN Generate nanowire Hamiltonian from bulk HR object
%
%   H_hr = HNANOWIRE_GEN(H_hr,Nslab,np,vacuum_mode,options) constructs
%   a nanowire Hamiltonian by repeating the unit cell along specified directions.
%
%   INPUT ARGUMENTS:
%       H_hr - Bulk HR object
%       Nslab - Supercell dimensions [Nx,Ny,Nz] (default: [1,1,1])
%       np - Number of parallel workers (default: 1)
%       vacuum_mode - Vacuum boundary mode (default: 0)
%       options - Structure with parameters:
%           symbolic: Use symbolic representation (logical, default: false)
%           fast: Use fast mode (logical, default: true)
%
%   OUTPUT ARGUMENTS:
%       H_hr - Nanowire HR object
%
%   NOTES:
%       - Supports parallel computation when np > 1
%       - Handles both sparse and dense matrix formats
%       - Automatically generates orbital positions
%
%   SEE ALSO:
%       HR, Hnum_list_wire_iz_gen
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]
%---------------------------
arguments
    H_hr HR;
    Nslab double= [1 1 1];
    np double{mustBeInteger} =1
    vacuum_mode double = 0;
    options.symbolic logical = false;
    options.fast logical = true;
end
if options.symbolic
    % not recommend
    H_hr.coe = true;
else
    H_hr.coe = false;
    H_hr.num = true;
end
%--------  reshape  --------
sort_dir = [3,2,1];
vectorList = double(H_hr.vectorL);
[vector_list_init,sort_label] = sortrows(vectorList,sort_dir) ;% sort fin_dir
H_hr = H_hr.reseq(':',sort_label);
WANNUM =  H_hr.WAN_NUM;
Ns = diag(Nslab);
H_hr.Rm = Ns * H_hr.Rm;
% --------  init  --------
WAN_NUM_x = WANNUM*Nslab(1);
WAN_NUM_y = WAN_NUM_x*Nslab(2);
% --------  3D begining  --------
[unique_z,unique_label_z]= unique(vector_list_init(:,3),'rows');
cutlist_z = HR.unique_label2cutlist(unique_label_z,H_hr.NRPTS);
NRPTS_z = length(unique_z);
vector_list_wire = (zeros(NRPTS_z,3));
vector_list_wire(:,3) = (unique_z);
% init
vertor_list_xy{NRPTS_z} = vector_list_init(cutlist_z(NRPTS_z,1):cutlist_z(NRPTS_z,2),:);
Hnum_list_wire{NRPTS_z} = sparse(WAN_NUM_y,WAN_NUM_y);
if strcmp(H_hr.Type,'sparse')
    Hnum_list_xy{NRPTS_z} = H_hr.HnumL(cutlist_z(NRPTS_z,1):cutlist_z(NRPTS_z,2));
else
    Hnum_list_xy{NRPTS_z} = H_hr.HnumL(:,:,cutlist_z(NRPTS_z,1):cutlist_z(NRPTS_z,2));
end
% --------  vector  --------
if  strcmp(H_hr.Type,'sparse')
    for iz = 1:NRPTS_z-1
        %                     vector_list_wire(iz,3) = unique_z(iz);
        vertor_list_xy{iz} = vector_list_init(cutlist_z(iz,1):cutlist_z(iz,2),:);
        Hnum_list_xy{iz} = H_hr.HnumL(cutlist_z(iz,1):cutlist_z(iz,2));
        Hnum_list_wire{iz} = sparse(WAN_NUM_y,WAN_NUM_y);
    end
else
    for iz = 1:NRPTS_z-1
        %                     vector_list_wire(iz,3) = unique_z(iz);
        vertor_list_xy{iz} = vector_list_init(cutlist_z(iz,1):cutlist_z(iz,2),:);
        Hnum_list_xy{iz}   = H_hr.HnumL(:,:,cutlist_z(iz,1):cutlist_z(iz,2));
        Hnum_list_wire{iz} = sparse(WAN_NUM_y,WAN_NUM_y);
    end
end

TYPE = H_hr.Type;
if np >1
    disp('parallel mode, we will use local settings, please set before.');
    np_handle = parpool('local',np);
    parfor iz = 1:NRPTS_z
        fprintf('Gen (%d/%d) NRPT z \n',iz,NRPTS_z);
        Hnum_list_xy_iz = Hnum_list_xy{iz};
        vertor_list_xy_iz = vertor_list_xy{iz};
        Hnum_list_wire{iz} = HR.Hnum_list_wire_iz_gen(Hnum_list_xy_iz,vertor_list_xy_iz,iz,WANNUM,WAN_NUM_x,WAN_NUM_y,Nslab,TYPE);
    end
    delete(np_handle);
else
    for iz = 1:NRPTS_z
        fprintf('Gen (%d/%d) NRPT z \n',iz,NRPTS_z);
        Hnum_list_xy_iz = Hnum_list_xy{iz};
        vertor_list_xy_iz = vertor_list_xy{iz};
        Hnum_list_wire{iz} = HR.Hnum_list_wire_iz_gen(Hnum_list_xy_iz,vertor_list_xy_iz,iz,WANNUM,WAN_NUM_x,WAN_NUM_y,Nslab, TYPE);
    end
end
if strcmp(H_hr.Type,'sparse')
    H_hr.HnumL = Hnum_list_wire;
else
    HnumL_temp = zeros(WAN_NUM_y,WAN_NUM_y,NRPTS_z);
    for iz = 1:NRPTS_z
        HnumL_temp(:,:,iz) = full(Hnum_list_wire{iz});
    end
    H_hr.HnumL = HnumL_temp ;
end
H_hr.vectorL = (vector_list_wire);
H_hr.orbL = H_hr.nanowire_orb(Nslab,vacuum_mode,'fast',options.fast);
for i = 1:H_hr.Dim
    if  Nslab(i)<1
        Nslab(i) = 1;
    end
end
num_sc = Nslab(1)* Nslab(2)* Nslab(3);
H_hr.elementL = repmat(H_hr.elementL,[num_sc,1]);
H_hr.quantumL = repmat(H_hr.quantumL,[num_sc,1]);
end