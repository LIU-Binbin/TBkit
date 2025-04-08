function HoutL = HCAR_gen(Hfun, klist_cart, Norb)
%HCAR_GEN Generate Hamiltonian matrices for given k-points
%   HoutL = HCAR_GEN(HFUN, KLIST_CART, NORB) computes Hamiltonian matrices
%   for each k-point in Cartesian coordinates
%
%   Inputs:
%   HFUN       - Function handle that returns Norb×Norb Hamiltonian matrix
%                for given k-point [kx, ky, kz]. Expected signature: 
%                @(kx,ky,kz) returning Norb×Norb matrix
%   KLIST_CART - N×3 array of k-point coordinates in Cartesian system
%   NORB       - Number of orbitals (determines matrix dimensions)
%
%   Output:
%   HOUTL      - 3D array of size Norb×Norb×N containing Hamiltonian matrices
%                for each k-point in KLIST_CART

    arguments
        Hfun function_handle;
        klist_cart double;
        Norb;
    end
    
    kn = size(klist_cart, 1);
    HoutL = zeros(Norb, Norb, kn);
    
    for i = 1:kn
        Input = num2cell(klist_cart(i, :));
        HoutL(:, :, i) = Hfun(Input{:});
    end
end