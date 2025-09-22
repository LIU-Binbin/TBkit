%% Example script: Constructing bismuth surface bands, building a moiré Hamiltonian,
%  and visualizing Fermi surfaces using TBkit.
%
%  This script demonstrates:
%    1. Building a tight-binding Hamiltonian with symmetry constraints.
%    2. Applying moiré potential terms.
%    3. Generating band structures and Fermi surfaces.
%    4. Constructing a moiré Hamiltonian and plotting its bands and Fermi surfaces.
clear
%% Define Pauli matrices in orbital (tau) and spin (sigma) space
sigma_0 = [1 0;0 1];
sigma_x = [0 1;1 0];
sigma_y = [0 -1i;1i 0];
sigma_z = [1 0;0 -1];
tau_0   = [1 0;0 1];
tau_x   = [0 1;1 0];
tau_y   = [0 -1i;1i 0];
tau_z   = [1 0;0 -1];

%% Initialize symmetry-constrained tight-binding model
%  4×4 model following the reference paper
Graphene = HR(4);
Graphene = Graphene<'POSCAR_spinful_0';

% Construct nearest-neighbor hoppings (from the paper)
Graphene = Graphene.nn([2,2,0],1e-15,10,"onsite",true);
sym(Graphene.nn_store);

% Cut the Hamiltonian up to 2nd order neighbors
Graphene2 = Graphene.init('level_cut',2,"onsite",true,'fast',true);

% Adjust lattice vectors for better precision
Graphene2.Rm = [1 0 0; -1/2 sqrt(3)/2 0; 0 0 10];

%% Generate basis functions and symmetry operators
Graphene.POSCAR_gen();
!mv POSCAR_gen.vasp POSCAR
BasisFunction = BasisFunc(Graphene2);

Tr  = Oper.time_reversal(3,diag([1 1 1 1]));
C3z = Oper.rotation(1/3,[0,0,1],false);
C3z = C3z.attachRm(Graphene2.Rm);
C3z.U = BasisFunction.rotation('Oper',C3z);

Mx  = Oper.mirror([1,0,0]);
Mx  = Mx.attachRm(Graphene2.Rm);
Mx.U = BasisFunction.rotation('Oper',Mx);

%% Define symmetry matrices for Bi(111)/FeCl2 system
Mx.U = blkdiag(-1i*sigma_z,[0 exp(1i*pi/3);exp(1i*2*pi/3) 0]); 
C3z.U = blkdiag(-sigma_0,[exp(-1i*pi/3) 0; 0 exp(1i*pi/3)]);
Tr.U  = blkdiag(1i*sigma_y,1i*sigma_y);

%% Apply symmetry constraints to the moiré potential
Vm = Graphene.nn([2,2,0],1e-20,200,"onsite",true);
Vm = Vm.init('level_cut',0,"onsite",true,'fast',true);
Vm = Vm.applyOper([C3z,Mx],'generator',true,'fast',true); % TR does not constrain Vm
Vmo = Vm.GenfromOrth('fromCvectorL',true);

%% Apply symmetry constraints to the Hamiltonian
Groups = generate_group([C3z,Mx,Tr]);
Graphene_test = Graphene2.applyOper([C3z,Mx,Tr], 'generator', true);
disp('Hamiltonian after symmetry constraints:');
list(Graphene_test);

Graphene_test2 = Graphene_test.GenfromOrth('fromCvectorL',true);
list(Graphene_test2);
Graphene_test2.show('HOPPING');

%% Band plot (parameter substitution)
gamma__r_1  = -0.5;  gamma__r_2  = 2.962;  gamma__r_3  = 0;     gamma__r_4  = -0.32;
gamma__r_5  = 0.658; gamma__r_6  = -0.658; gamma__r_7  = -0.32;
gamma__r_8  = 0;     gamma__r_9  = 2.962;  gamma__r_10 = 0;     gamma__r_11 = 0; 
gamma__r_12 = -1.2;

Graphene_test2_n = Graphene_test2.Subsall();
Graphene_test2_n = Graphene_test2_n<'KPOINTS_hex_2D';
Graphene_test2_n.bandplot([-1,3]);

%% Symmetry check in symbolic form
syms k_x k_y k_z real;
simplify(subs(Tr.U*conj(Graphene_test2_n.sym)/(Tr.U),[k_x k_y k_z],[-k_x -k_y -k_z])-Graphene_test2_n.sym)
S1 = simplify(subs(C3z.U*(Graphene_test2_n.sym)/(C3z.U),[k_x k_y k_z],[-k_x/2-sqrt(3)*k_y/2 -k_y/2+sqrt(3)*k_x/2 k_z])-Graphene_test2_n.sym);
S2 = simplify(subs(Mx.U*(Graphene_test2_n.sym)/(Mx.U),[k_x k_y k_z],[-k_x k_y k_z])-Graphene_test2_n.sym);
round(subs([S1],[k_x k_y k_z],[1.23 2.45 1.93])) 
round(subs([S2],[k_x k_y k_z],[0.13 3.26 2.23]))

%% Generate 3D Fermi surface for the original band
kmesh = [151 151];
lk = 2;
[EIGENCAR_3D,klist1,klist2] = Graphene_test2_n.EIGENCAR_gen_3D(kmesh,[[-lk -lk 0]/2; lk*[1 0 0]; lk*[0 1 0]]);
R = 0.25; % mBZ / BZ ratio

Efermi = 0.0;
[~,Ax]= Figs(1);
bandplot3d(EIGENCAR_3D-Efermi,klist1,klist2,'Ecut',[0 0.01],'ax',Ax,'xlabel','k_x','ylabel','k_y','title',Efermi);

%% Construct moiré Hamiltonian
Graphene_test2_n = Graphene_test2_n<'KPOINTS_hex_2D_MmKmGKmMm'; 
R = 0.25; % mBZ/BZ ratio
N = Graphene_test2_n.WAN_NUM;
klist_cart = Graphene_test2_n.klist_cart * R;

% Define moiré coupling potential
v0 = 0.012;
Vm = subs(Vmo.sym, Vmo.symvar_list, [v0, v0, v0]);

% Construct moiré Hamiltonian function
Hm = construct_moire_H(Graphene_test2_n,Vm,0.25);

% Diagonalize Hm across k-mesh
numk = size(klist_cart,1);
EIGENCAR = zeros(N*7, numk);
clear WAVECAR EIGENCAR
for j = 1:numk
    kpt = klist_cart(j,:);
    Hmk = Hm(kpt(1),kpt(2),kpt(3));
    [A, U] = eig(Hmk);
    [A, U] = park.sorteig(U, A);
    WAVECAR(:,:,j) = A;
    EIGENCAR(:, j) = diag(U);
end

% Band plot for moiré bands
klist_l     = Graphene_test2_n.klist_l;
kpoints_l   = Graphene_test2_n.kpoints_l;
kpoints_name= Graphene_test2_n.kpoints_name;
bandplot(EIGENCAR, [-0.2, 0.1], klist_l, kpoints_l, kpoints_name);

%% Fermi surface calculation for moiré bands
kmesh = [151, 151];
method = 'area';
radius = 0.0129;
lk    = 3;
blk = 7; %num of diagoanl block 

if method == 'area'
    [klist_cart, klist_frac, klist_r_plot, sizemesh, Gk_, Grid] = ...
        vasplib.kmesh2D(Graphene_test2_n.Rm/R, ...
        'knum1', kmesh(1), 'knum2', kmesh(2), ...
        'kstart', [-lk, -lk, 0]/2, ...
        'kdir1', lk*[1 0 0], ...
        'kdir2', lk*[0 1 0]);
end

clear WAVECAR EIGENCAR
numk = size(klist_cart,1);
EIGENCAR = zeros(blk*N, numk);
Hm = construct_moire_H(Graphene_test2_n,Vm,0.25);
for j = 1:numk
    kpt = klist_cart(j,:);
    Hmk = Hm(kpt(1),kpt(2),kpt(3));
    Hmk = (Hmk + Hmk')/2;
    [A, U] = eig(Hmk);
    [A, U] = park.sorteig(U, A);
    WAVECAR(:,:,j)=A;
    EIGENCAR(:, j)=diag(U);
end

% Reshape eigenvalues into 3D array [kx, ky, band]
EIGENCAR_3Dsh = reshape(EIGENCAR.', kmesh(1), kmesh(2), []);

%% Plot moiré Fermi surface
Efermi = 0.0;
bandplot3d(EIGENCAR_3Dsh-Efermi, ...
    reshape(klist_cart(:,1),kmesh), ...
    reshape(klist_cart(:,2),kmesh), ...
    'Ecut',[-0.0020,0.002], 'xlabel',Efermi);

%% construct_moire_H.m
% Construct the Moiré Hamiltonian from a tight-binding model from HR class.
%
% This function builds the Moiré Hamiltonian Hm given a parent Hamiltonian HR,
% a Moiré potential matrix Vmijs, and mBZ/BZ ratio R.
% The resulting Hamiltonian includes the central block and six surrounding
% translated copies (±Gm shifts), along with coupling terms that connect them, if restore_C3z = True.
%
% Usage:
%   Hm = construct_moire_H(HR, Vmijs, R, restore_C3z)
%
% Inputs:
%   HR           - HR object from TBkit, containing the parent TB Hamiltonian
%   Vmijs        - Moiré coupling potential matrix of size (N × N), typically
%                  block-diagonal, e.g. blkdiag(v1*I, v2*I).
%   R            - Scaling factor between mBZ and BZ (default: 0.25).
%   restore_C3z  - Logical flag, if true, restore C3z rotational symmetry by
%                  coupling between neighboring satellite blocks (default: true).
%
% Outputs:
%   Hm           - Function handle @(kx,ky,kz) returning the Moiré Hamiltonian
%                  matrix of dimension N*(1+6), i.e. center + 6 neighbors.
%
% Notes:
%   - The function builds symbolic expressions in k-space, then converts them
%     into a MATLAB function handle using matlabFunction.
%   - Coupling terms are inserted between central and satellite Hamiltonians,
%     and optionally between satellites to enforce C3 symmetry.
%   - Hermiticity is enforced explicitly at the end.

function Hm = construct_moire_H(HR, Vmijs, R, restore_C3z)
arguments
    HR
    Vmijs
    R                  % Ratio mBZ/BZ
    restore_C3z logical = true;
end

% Extract reciprocal lattice vectors from parent HR
Gk = HR.Gk;

% Construct Moiré reciprocal lattice vectors
Gm = R * Gk; 
Gm(3,:) = -Gm(1,:) + Gm(2,:);   % Third vector for hexagonal symmetry
N = HR.WAN_NUM;                 % Number of Wannier orbitals per unit cell

% Define symbolic momentum variables
syms k_x k_y k_z real;

% Parent Hamiltonian as function of k
Hfun = matlabFunction(HR.sym, 'Vars', [k_x, k_y, k_z]);

% Extend Gm to include both +Gm and -Gm shifts
Gm = [Gm; -Gm];

% Construct block-diagonal Hamiltonian: center + 6 neighbors
Hblocks = cell(7,1);
Hblocks{1} = Hfun(k_x, k_y, k_z);  % central block
for i = 1:6
    k_shift = [k_x, k_y, k_z] - Gm(i,:);
    Hblocks{i+1} = Hfun(k_shift(1), k_shift(2), k_shift(3));
end
Hm = blkdiag(Hblocks{:}); 

% Add couplings between central block and each satellite
for i = 1:6
    row_idx = 1:N;
    col_idx = (N*i + 1):(N*(i+1));
    
    % Insert Hermitian coupling
    Hm(row_idx, col_idx) = Vmijs;
    % Hm(col_idx, row_idx) = Vmijs';
end

% Optionally restore C3z symmetry by adding couplings between satellites
for i = 1:6
    j = mod(i,6) + 1;  % next neighbor index (1..6 cyclic)
    idx_i = (N*i + 1):(N*(i+1));
    idx_j = (N*j + 1):(N*(j+1));
    
    if restore_C3z
        Hm(idx_i, idx_j) = Vmijs;
        Hm(idx_j, idx_i) = Vmijs';
    end
end

% Enforce Hermiticity explicitly
Hm = (Hm + Hm')/2;

% Convert symbolic Hamiltonian into function handle @(kx,ky,kz)
Hm = matlabFunction(Hm, 'Vars', [k_x, k_y, k_z]);
end
