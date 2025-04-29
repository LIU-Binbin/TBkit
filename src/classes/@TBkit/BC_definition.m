function BC = BC_definition(Hsym,para1,para2,epsilon,options)
%BC_DEFINITION Calculate Berry curvature from symbolic Hamiltonian
%
%   Syntax:
%       BC = BC_definition(Hsym,para1,para2,epsilon,options)
%
%   Description:
%       Computes Berry curvature using eigenvector derivatives for a
%       symbolic Hamiltonian. Based on the definition:
%       Ω = ∂A₂/∂para1 - ∂A₁/∂para2 where A is the Berry connection.
%
%   Inputs:
%       Hsym      - Symbolic Hamiltonian matrix
%       para1/2   - Parameters to differentiate with respect to
%       epsilon   - Scaling factor (default=1)
%       options   - Options structure:
%                   BAND_index - Band indices to include
%
%   Output:
%       BC - Berry curvature values for selected bands
arguments
    Hsym sym;
    para1 sym;
    para2 sym;
    epsilon = 1;
    options.BAND_index = [];
end
[Norb,~] = size(Hsym);
set_divide = 2;
if isempty(options.BAND_index)
    BAND_index = 1:(Norb/set_divide);
else
    BAND_index = options.BAND_index;
end
%HsymP1 = diff(Hsym,para1);
%HsymP2 = diff(Hsym,para2);
[W,~] = eig(Hsym);
W = W(:,BAND_index);
[~,Nbands] = size(W);
BC = zeros([Nbands,1],class(W));
%u = zeros(size(W),class(W));
for i = 1:Nbands
    Eigenvector_sym = (TBkit.NomalizeEigenvector(W(:,i)));
    A_1 = TBkit.BerryConnection_definition(Eigenvector_sym,para1);
    A_2 = TBkit.BerryConnection_definition(Eigenvector_sym,para2);
    BC(i) = TBkit.BerryCurvature_definition(A_1,A_2,para1,para2);
    %WL(:,:,i) = HsymP1*u(:,i)*u(:,i)'*HsymP2;
end
BC = BC*epsilon;
%E = simplify(diag(U));
% $\begin{aligned} \Omega_{\mu \nu}^{n}(\boldsymbol{R}) &=\frac{\partial}{\partial R_{\mu}} A_{\nu}^{n}(\boldsymbol{R})-\frac{\partial}{\partial R_{\nu}} A_{\mu}^{n}(\boldsymbol{R}) \\ &=i\left[\frac{\partial}{\partial R_{\mu}}\left\langle n(\boldsymbol{R})\left|\frac{\partial}{\partial R_{\nu}}\right| n(\boldsymbol{R})\right\rangle-(\nu \leftrightarrow \mu)\right] \end{aligned}$
%Bc = simplify(diff(A_2,para1)-diff(A_1,para2));
end