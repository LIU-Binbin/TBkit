function [EIGENCAR,orb_list,WAVECAR] = EIGENCAR_gen_disk(H_hr,Nslab,fermi,norb_enforce,kpoints,vacuum_mode,np)
% EIGENCAR_GEN_DISK Generate eigenstate data for disk/nanowire geometry
%
%   [EIGENCAR,orb_list,WAVECAR] = EIGENCAR_GEN_DISK(H_hr,Nslab,fermi,norb_enforce,kpoints,vacuum_mode,np)
%   calculates eigenstates for a disk or nanowire geometry derived from the bulk Hamiltonian.
%
%   INPUT ARGUMENTS:
%       H_hr - Bulk Hamiltonian in HR format
%       Nslab - Nanowire dimensions [Nx,Ny,Nz] (default: [10,10,0])
%       fermi - Fermi level (default: 0)
%       norb_enforce - Number of bands to enforce (-1 for all, default: -1)
%       kpoints - k-points to calculate (default: [0,0,0])
%       vacuum_mode - Vacuum boundary mode (default: 1)
%       np - Parallel processing flag (default: 0)
%
%   OUTPUT ARGUMENTS:
%       EIGENCAR - Eigenvalue data
%       orb_list - Orbital list
%       WAVECAR - Wavefunction data
%
%   NOTES:
%       - Can handle parallel processing when np > 1
%       - Automatically generates nanowire Hamiltonian if Nslab specified
%
%   SEE ALSO:
%       HR, Hnanowire_gen
%
%   AUTHOR:
%       [Your Name] ([Your Email])
%       [Creation Date]

import TBkit_tool.*
if nargin < 2
Nslab = [10 10 0];
end
if nargin < 3
fermi = 0;
end
if nargin < 4
norb_enforce  = -1;
end
if nargin < 5
kpoints = [0 0 0];
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
H_hr_wire = H_hr_wire.sparse();
Hnum_list_wire = H_hr_wire.HnumL;
vector_list_wire = double(H_hr_wire.vectorL);
[nz,~] = size(vector_list_wire);
factor_list_wire = exp(1i*2*pi*kpoints*vector_list_wire');
orb_list  = H_hr_wire.orbL;
NWAVE = H_hr_wire.WAN_NUM;
if norb_enforce <0
NBANDS=NWAVE;
elseif norb_enforce >0
NBANDS=norb_enforce;
else
end
kn = size(kpoints,1);
EIGENCAR = zeros(NBANDS,kn);
WAVECAR = zeros(NWAVE,NBANDS);
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
[A, U]= HR.sorteig(U,A);
else
end
EIGENCAR(:,ki) = diag(U);
WAVECAR= A;
fprintf('%d th kpoint(%7.4f, %7.4f, %7.4f) has been calculated in %d kpoints total\n',...
ki,kpoints(1),kpoints(2),kpoints(3),kn);
end
if np >1
delete(np_handle);
end
end
