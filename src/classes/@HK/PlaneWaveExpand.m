function H_hk_out = PlaneWaveExpand(H_hk, N, ExpandDirection)
%PLANEWAVEEXPAND Expand Hamiltonian in plane wave basis
%
% Syntax:
%   H_expanded = PlaneWaveExpand(H_hk)
%   H_expanded = PlaneWaveExpand(H_hk, N)
%   H_expanded = PlaneWaveExpand(H_hk, N, ExpandDirection)
%
% Inputs:
%   H_hk - Original HK Hamiltonian
%   N - Number of plane waves per direction (default=5)
%   ExpandDirection - Directions to expand [x y z] (default=[1 1 0])
%
% Output:
%   H_hk_out - Expanded Hamiltonian with:
%     - Basis_num increased by (2N+1)^d (d=expansion dimensions)
%     - Modified reciprocal space terms
%
% Description:
%   Expands the Hamiltonian by introducing plane wave components:
%   1. Creates reciprocal lattice vectors up to order N
%   2. Modifies k-space terms with G-vector shifts
%   3. Handles both symbolic and numeric coefficient cases
%
% Note:
%   Automatically converts numeric coefficients to symbolic
%   Expansion preserves original Hamiltonian structure
%
% Example:
%   H_exp = Hk.PlaneWaveExpand(3,[1 1 1]); % 3D expansion up to ±3G
arguments
    H_hk HK;
    N double{mustBeInteger} = 5;
    ExpandDirection = [1 1 0];
end
if H_hk.num
    H_hk.HcoeL = sym(H_hk.HnumL);
    H_hk.num = false; H_hk.coe = true;
else
end
ExpandNum =(2*N+1)^sum(ExpandDirection) ;
BasisNum = ExpandNum * H_hk.Basis_num;
NL = (-N:N).';
if ExpandDirection(1)
    mL = NL;
else
    mL = 0;
end
if ExpandDirection(2)
    nL = NL;
else
    nL = 0;
end
if ExpandDirection(3)
    lL = NL;
else
    lL = 0;
end
mnlL = [...
    kron(mL,ones(length(nL)*length(lL),1)),...
    kron(ones(length(mL),1),kron(nL,ones(length(lL),1))),...
    kron(ones(length(mL)*length(nL),1),lL)];
pqoL = mnlL;
syms G_0_0_0 k_x k_y k_z real;
H_hk.HcoeL = H_hk.HcoeL;
HcoeListpre = H_hk.HcoeL;
pqoL_G = kron(ones(H_hk.Basis_num ,1 ),pqoL * H_hk.Gk);
Symvarlist = symvar(HcoeListpre);
SymvarlistStrL = string(Symvarlist);
G_Symvarlist = Symvarlist(contains(SymvarlistStrL,"G_"));
nG_Symvarlist = length(G_Symvarlist);
G_SymvarCell{nG_Symvarlist} = sym('0','real');
G_SymvarCellSubs{nG_Symvarlist} = zeros(ExpandNum);
for i = 1:nG_Symvarlist
    G_SymvarCell{i} = G_Symvarlist(i);
    VectorTmp = park.Variable2Vector(G_Symvarlist(i));
    G_SymvarCellSubs{i} = park.DeltaList(mnlL,VectorTmp+pqoL);
end
G_SymvarCell{nG_Symvarlist+1} = G_0_0_0;
G_SymvarCellSubs{nG_Symvarlist+1} = eye(ExpandNum);
HcoeListpre = expand(simplify(HcoeListpre*G_0_0_0));
HcoeListpreStr = string(HcoeListpre);
HcoeListpreStr = strrep(HcoeListpreStr,"G_0_0_0*G","G");
HcoeListpreModify = str2sym(HcoeListpreStr);
HcoeList = subs(HcoeListpreModify,G_SymvarCell,G_SymvarCellSubs);
H_hk_out = H_hk;
H_hk_out.HcoeL = HcoeList;
H_hk_out.Basis_num = BasisNum;
H_hk_symL = H_hk.HsymL;
H_hk_symL_2 = sym([ones(BasisNum,H_hk.Kinds)]);
for i = 2:H_hk.Kinds
    for j = 1:BasisNum
        H_hk_symL_2(j,i) = subs(H_hk_symL(i),[k_x,k_y,k_z],[...
            k_x+pqoL_G(j,1),...
            k_y+pqoL_G(j,2),...
            k_z+pqoL_G(j,3)...
            ]);
    end
end
H_hk_out.HsymL = H_hk_symL_2;
end
