function [Gwl,Gwb,Gwr] = GW_iter(H00,H01,w,eta,mu_max,infinity_small)
%GW_ITER Compute Green's functions using recursive algorithm
%
%   Syntax:
%       [Gwl,Gwb,Gwr] = GW_iter(H00,H01,w,eta,mu_max,infinity_small)
%
%   Description:
%       Implements the Sancho-Rubio recursive Green's function method
%       for surface and bulk Green's functions in layered systems.
%
%   Inputs:
%       H00           - Intra-layer Hamiltonian
%       H01           - Inter-layer coupling
%       w             - Frequency/energy point
%       eta           - Broadening parameter
%       mu_max        - Maximum recursion steps (default=100)
%       infinity_small - Convergence threshold (default=1e-10)
%
%   Outputs:
%       Gwl - Left surface Green's function
%       Gwb - Bulk Green's function
%       Gwr - Right surface Green's function
%
%   References:
%       Sancho et al., J. Phys. F: Met. Phys. 15 (1985) 851-858
%
%   See also: GREENCAR_gen, Tmatrix_iter
%%
if nargin <6
    infinity_small= 1e-10;
end
%%
if nargin <5
    mu_max = 100;
end
epsilon0 = H00;
alpha0 = H01;
beta0 = H01';

%oumega_bk = E_w;
oumega =  (w+1i*eta)*eye(length(H00));
%save('Green_prepare.mat');
epsilon_s0 = epsilon0;
epsilon_s_bar0 = epsilon0;
%[epsilon,epsilon_s,epsilon_s_bar,~,~] = GW_iter_coff(mu,oumega,alpha0,beta0,epsilon0);

for i = 1:mu_max
    Ml = (oumega - epsilon0)\alpha0;
    Mr = (oumega - epsilon0)\beta0;
    M1 = alpha0 * (Ml);
    M2 = beta0  * (Mr);
    M3 = alpha0 * (Mr);
    M4 = beta0  * (Ml);

    difference = M3+ M4;
    if norm(difference,'fro') < infinity_small
        %disp('reach_accuracy')
        epsilon0 = epsilon0 +M3+ M4;
        epsilon_s0 = epsilon_s0 + M3;
        epsilon_s_bar0 = epsilon_s_bar0 + M4;
        break;
    end
    %     disp(max( difference(:)));
    alpha0 = M1;
    beta0 = M2;
    epsilon0 = epsilon0 +M3+ M4;
    epsilon_s0 = epsilon_s0 + M3;
    epsilon_s_bar0 = epsilon_s_bar0 + M4;

    %     [epsilon,epsilon_s,epsilon_s_bar,alpha,beta,difference] = GW_iter_coff_once(oumega,epsilon0,epsilon_s0,epsilon_s_bar0,alpha0,beta0);
    %     alpha0 = alpha;
    %     beta0 = beta;
    %     epsilon0 = epsilon;
    %     epsilon_s0 = epsilon_s;
    %     epsilon_s_bar0 = epsilon_s_bar;
    %disp(i);
    %disp(max(difference(:)));

end
% epsilon = epsilon0;
% epsilon_s = epsilon_s0 ;
% epsilon_s_bar = epsilon_s_bar0;
Gwl = inv((oumega-epsilon_s0));
Gwb = inv((oumega-epsilon0));
Gwr = inv((oumega-epsilon_s_bar0));
end