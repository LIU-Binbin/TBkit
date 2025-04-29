function F = BerryPhaseLine_definition_sym_i(Hsym,options)
%BERRYPHASELINE_DEFINITION_SYM_I Calculate Berry phase line integral symbolically
%
%   Syntax:
%       F = BerryPhaseLine_definition_sym_i(Hsym,options)
%
%   Description:
%       Computes the Berry phase line integral for a symbolic Hamiltonian
%       using the connection 1-form A = i<ψ|∂ψ/∂k>·dk along a path.
%
%   Inputs:
%       Hsym      - Symbolic Hamiltonian matrix
%       options   - Options structure:
%                   BAND_index - Band indices to include (default=all)
%
%   Output:
%       F - Berry phase line integral values for selected bands
arguments
    Hsym sym;
    options.BAND_index = [];
end
% temp 2d
syms k_x k_y k_z delta_k_x delta_k_y delta_k_z real;
[Norb,~] = size(Hsym);
set_divide = 2;
if isempty(options.BAND_index)
    BAND_index = 1:(Norb/set_divide);
else
    BAND_index = options.BAND_index;
end
[W,~] = eig(Hsym);
W = W(:,BAND_index);
[~,Nbands] = size(W);
F = zeros([Nbands,1],class(W));
for i = 1:Nbands
    Eigenvector_sym = (TBkit.NomalizeEigenvector(W(:,i)));
    A_x = TBkit.BerryConnection_definition(Eigenvector_sym,k_x);
    A_y = TBkit.BerryConnection_definition(Eigenvector_sym,k_y);
    A_z = TBkit.BerryConnection_definition(Eigenvector_sym,k_z);
    F(i) = ([delta_k_x delta_k_y delta_k_z]*[A_x,A_y,A_x].');
    %WL(:,:,i) = HsymP1*u(:,i)*u(:,i)'*HsymP2;
end
end