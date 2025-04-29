function Bc = BC_kubo_sym(Hsym,para1,para2,epsilon,options)
%BC_KUBO_SYM Calculate Berry curvature using Kubo formula
%
%   Syntax:
%       Bc = BC_kubo_sym(Hsym,para1,para2,epsilon,options)
%
%   Description:
%       Computes Berry curvature using the Kubo formula approach for
%       symbolic Hamiltonians. Implements the formula:
%       Ω = iΣ_{n≠m} (<n|∂H/∂para1|m><m|∂H/∂para2|n> - (1↔2))/(E_n-E_m)^2
%
%   Inputs:
%       Hsym      - Symbolic Hamiltonian matrix
%       para1/2   - Parameters to differentiate with respect to
%       epsilon   - Scaling factor (default=1)
%       options   - Options structure:
%                   BAND_index - Band indices to include
%
%   Output:
%       Bc - Berry curvature values for selected bands
if nargin < 4
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
HsymP1 = diff(Hsym,para1);
HsymP2 = diff(Hsym,para2);
[W,U] = eig(Hsym);
E = simplify(diag(U));
WL = zeros(Norb,Norb,Norb,class(W));
%WL2 = WL;% full
for i = 1:Norb
    W(:,i) = (TBkit.NomalizeEigenvector(W(:,i)));
    WL(:,:,i) = HsymP1*W(:,i)*W(:,i)'*HsymP2;
    %WL2(:,:,i) = HsymP2*W(:,i)*W(:,i)'*HsymP1;% full
end
u = W(:,BAND_index);
[Norb,Nbands] = size(u);
Bc = zeros([Nbands,1],class(W))+1i*zeros([Nbands,1],class(W));
for i = 1:Nbands
    for j = 1:Norb
        %$\Omega_{\mu \nu}^{n}(\boldsymbol{R})=i \sum_{n^{\prime} \neq n} \frac{\left\langle n\left|\frac{\partial H}{\partial R_{\mu}}\right| n^{\prime}\right\rangle\left\langle n^{\prime}\left|\frac{\partial H}{\partial R_{\nu}}\right| n\right\rangle-(\nu \leftrightarrow \mu)}{\left(E_{n}-E_{n^{\prime}}\right)^{2}}$
        if j ~= i && simplify(E(i)-E(j)) ~= sym(0)
            %                         BC(i) = BC(i) + (1i/(E(i)-E(j)))*((u(:,i)'*HsymP1*u(:,i))*(u(:,i)'*HsymP2*u(:,i))- ...
            %                 (u(:,i)'*HsymP2*u(:,j))*(u(:,j)'*HsymP1*u(:,i)));
            % if Hermite
            Bc(i) = Bc(i) +(1/(E(i)-E(j))^2)*(u(:,i)'*WL(:,:,j)*u(:,i));
            % full
            %Bc(i) = Bc(i)+(1/(E(i)-E(j))^2)*((u(:,i)'*WL(:,:,j)*u(:,i)) -(u(:,i)'*WL2(:,:,j)*u(:,i)) ) ;
        end
    end
end
Bc = (Bc*2i*epsilon);
%Bc = (Bc*1i*epsilon);% full
end