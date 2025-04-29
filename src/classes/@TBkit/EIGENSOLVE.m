function [EIGENCAR,WAVECAR,HoutL] = EIGENSOLVE(Hfun,klist_cart,Norb,opt)
%EIGENSOLVE Solve eigenvalue problem for Hamiltonian across k-points
%
%   Syntax:
%       [EIGENCAR,WAVECAR,HoutL] = EIGENSOLVE(Hfun,klist_cart,Norb,opt)
%
%   Description:
%       Diagonalizes Hamiltonian at multiple k-points to compute eigenvalues
%       and eigenvectors. Supports both function handles and TBkit objects.
%
%   Inputs:
%       Hfun       - Hamiltonian function handle or TBkit object
%       klist_cart - Array of k-points in Cartesian coordinates
%       Norb       - Number of orbitals (default=from Hfun)
%       opt        - Options structure:
%                    Hermitian - Enforce Hermiticity (default=true)
%
%   Outputs:
%       EIGENCAR - Eigenvalues (bands × k-points)
%       WAVECAR  - Wavefunctions (orbitals × bands × k-points)
%       HoutL    - Hamiltonian matrices (orbitals × orbitals × k-points)
%
%   See also: HCAR_gen, TBkit
arguments
    Hfun
    klist_cart = [];
    Norb = -1;
    opt.Hermitian = true;
end
if isa(Hfun,'TBkit')
    if ~isempty(Hfun.klist_cart)
        klist_cart = Hfun.klist_cart;
    end
    Norb = Hfun.Basis_num;
    Hfun = Hfun.Hfun;
end
if Norb < 0
    Norb = length(Hfun(0,0,0));
end
%HoutL = TBkit.HCAR_gen(Hfun,klist_cart,Norb);
if isempty(klist_cart)
    warning('Empty klist!');
    EIGENCAR = [];
    WAVECAR = [];
    HoutL = [];
    return;
end
kn = size(klist_cart,1);
EIGENCAR = zeros(Norb,kn);
if nargout >1
    WAVECAR  = zeros(Norb,Norb,kn);
end
if nargout >2
    HoutL = zeros(Norb,Norb,kn);
end
for i = 1:kn
    Input = num2cell(klist_cart(i,:));
    Hout= Hfun(Input{:});
    if opt.Hermitian
        Hout = (Hout+Hout')/2;
    end
    [A, U]=eig(Hout);
    [A, U]= park.sorteig(U,A);
    EIGENCAR(:,i) = diag(U);
    if nargout >1
        WAVECAR(:,:,i)=A;
    end
    if nargout >2
        HoutL(:,:,i) = Hout;
    end
end

end