function [EIGENCAR,orb_list,WEIGHTCAR,klist_l,kpoints_l,kpoints_name] = EIGENCAR_gen_wire(H_hr,Nslab,fermi,norb_enforce,KPOINTS_wire,vacuum_mode,np)

% EIGENCAR_GEN_WIRE Generate eigenstate data for nanowire geometry
%
%   [EIGENCAR,orb_list,WEIGHTCAR,klist_l,kpoints_l,kpoints_name] = EIGENCAR_GEN_WIRE(H_hr,Nslab,fermi,norb_enforce,KPOINTS_wire,vacuum_mode,np)
%   calculates eigenstates for a nanowire derived from bulk Hamiltonian.
%
%   INPUT ARGUMENTS:
%       H_hr - Bulk Hamiltonian in HR format
%       Nslab - Nanowire dimensions [Nx,Ny,Nz] (default: [0,0,0] for bulk)
%       fermi - Fermi level (default: 0)
%       norb_enforce - Number of bands to enforce (-1 for all, default: -1)
%       KPOINTS_wire - k-point path file or string (default: 'KPOINTS_wire')
%       vacuum_mode - Vacuum boundary mode (default: 1)
%       np - Number of parallel workers (default: 0, no parallel)
%
%   OUTPUT ARGUMENTS:
%       EIGENCAR - Eigenvalue data (NBANDS x kn array)
%       orb_list - Orbital list
%       WEIGHTCAR - Weight data for hinge states
%       klist_l - k-point list in lattice coordinates
%       kpoints_l - k-point coordinates
%       kpoints_name - k-point labels
%
%   NOTES:
%       - Supports parallel computation when np > 1
%       - Automatically generates nanowire Hamiltonian if Nslab specified
%       - Uses sparse matrix operations for efficiency
%
%   SEE ALSO:
%       HR, Hnanowire_gen, COLORCAR_gen, sorteig
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

import TBkit_tool.*
if nargin < 2
Nslab = [0 0 0];
end
if nargin < 3
fermi = 0;
end
if nargin < 4
norb_enforce  = -1;
end
if nargin < 5
KPOINTS_wire = 'KPOINTS_wire';
end
if nargin < 6
vacuum_mode = 1;
end
if nargin < 7
np = 0;
end
if isequal(Nslab, [0,0,0])
H_hr_wire = H_hr;
else
H_hr_wire = H_hr.Hnanowire_gen(Nslab,np,vacuum_mode);
end
H_hr_wire = H_hr_wire < KPOINTS_wire;
H_hr_wire = H_hr_wire.sparse();
if np >1
disp('parallel mode, we will use local settings, please set before.');
np_handle = parpool('local',np);
else
end
Hnum_list_wire = H_hr_wire.HnumL;
vector_list_wire = double(H_hr_wire.vectorL);
[nz,~] = size(vector_list_wire);
factor_list_wire = exp(1i*2*pi*H_hr_wire.klist_frac*vector_list_wire');
orb_list  = H_hr_wire.orbL;
HSVCAR_hinge = TBkit.HSVCAR_gen(orb_list,'hinge');
NWAVE = H_hr_wire.WAN_NUM;
if norb_enforce <0
NBANDS=NWAVE;
elseif norb_enforce >0
NBANDS=norb_enforce;
else
end
kn = size(H_hr_wire.klist_frac,1);
EIGENCAR = zeros(NBANDS,kn);
WEIGHTCAR = zeros(NBANDS,kn);
if np >1 && kn >1
fprintf('begining parfor loop, we have %d kpoints\n',kn);
parfor ki =1:kn
A = zeros(NWAVE,NWAVE);
U = zeros(NWAVE,NWAVE);
factor_list = factor_list_wire(ki,:);
Hout = sparse(NWAVE,NWAVE);
for iz = 1:nz
Hout = Hout+Hnum_list_wire{iz}*factor_list(iz);
end
Hout = (Hout+Hout')/2;
if norb_enforce <0
[A, U]=eig(full(Hout));
elseif norb_enforce >0
[A, U]=eigs(Hout,NBANDS,fermi);
[A, U]= park.sorteig(U,A);
else
end
EIGENCAR(:,ki) = diag(U);
[~,WEIGHTCAR(:,ki)] = COLORCAR_gen(A,HSVCAR_hinge);
fprintf('%d th kpoints has been calculated in %d kpoints total\n',ki,kn);
end
elseif  kn <2
WEIGHTCAR = zeros(NWAVE,NBANDS);
for ki =1:kn
factor_list = factor_list_wire(ki,:);
Hout = sparse(NWAVE,NWAVE);
for iz = 1:nz
Hout = Hout+Hnum_list_wire{iz}*factor_list(iz);
end
Hout = (Hout+Hout')/2;
if norb_enforce <0
[A, U]=eig(full(Hout));
elseif norb_enforce >0
[A, U]=eigs(Hout,NBANDS,fermi);
[A, U]= park.sorteig(U,A);
else
end
EIGENCAR(:,ki) = diag(U);
if kn >1
[~,WEIGHTCAR(:,ki)] = COLORCAR_gen(A,HSVCAR_hinge);
else
WEIGHTCAR= A;
end
fprintf('%d th kpoint(%7.4f, %7.4f, %7.4f) has been calculated in %d kpoints total\n',...
ki,H_hr_wire.klist_frac(ki,1),H_hr_wire.klist_frac(ki,2),H_hr_wire.klist_frac(ki,3),kn);
end
else
for ki =1:kn
factor_list = factor_list_wire(ki,:);
Hout = sparse(NWAVE,NWAVE);
for iz = 1:nz
Hout = Hout+Hnum_list_wire{iz}*factor_list(iz);
end
Hout = (Hout+Hout')/2;
if norb_enforce <0
[A, U]=eig(full(Hout));
elseif norb_enforce >0
[A, U]=eigs(Hout,NBANDS,fermi);
[A, U]= park.sorteig(U,A);
else
end
EIGENCAR(:,ki) = diag(U);
[~,WEIGHTCAR(:,ki)] = COLORCAR_gen(A,HSVCAR_hinge);
fprintf('%d th kpoint(%7.4f, %7.4f, %7.4f) has been calculated in %d kpoints total\n',...
ki,H_hr_wire.klist_frac(ki,1),H_hr_wire.klist_frac(ki,2),H_hr_wire.klist_frac(ki,3),kn);
end
end
if np >1
delete(np_handle);
end
[klist_l,kpoints_l,kpoints_name] = H_hr_wire.kpath_information();
end
