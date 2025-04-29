function [Gwl,Gwb,Gwr] = Tmatrix_iter(H00,H01,w,eta,mu_max,infinity_small)
% TMATRIX_ITER Calculate Green's functions using T-matrix iteration method
%
% Syntax:
%   [Gwl,Gwb,Gwr] = Tmatrix_iter(H00,H01,w,eta)
%   [Gwl,Gwb,Gwr] = Tmatrix_iter(H00,H01,w,eta,mu_max)
%   [Gwl,Gwb,Gwr] = Tmatrix_iter(H00,H01,w,eta,mu_max,infinity_small)
%
% Description:
%   Computes left (Gwl), bulk (Gwb), and right (Gwr) Green's functions using
%   an iterative T-matrix approach for a given Hamiltonian.
%
% Input Arguments:
%   H00 - On-site Hamiltonian matrix
%   H01 - Nearest neighbor hopping matrix
%   w - Energy value (can be complex)
%   eta - Imaginary part for regularization (broadening)
%   mu_max - Maximum number of iterations (default: 100)
%   infinity_small - Convergence threshold (default: 1e-10)
%
% Output Arguments:
%   Gwl - Left surface Green's function
%   Gwb - Bulk Green's function
%   Gwr - Right surface Green's function
%
% Example:
%   H00 = diag([1,1]); H01 = eye(2);
%   [Gwl,Gwb,Gwr] = Tmatrix_iter(H00,H01,0.5,0.01);
if nargin <6
    infinity_small= 1e-10;
end
%%
if nargin <5
    mu_max = 100;
end
oumega =  (w+1i*eta)*eye(length(H00));
oumega_00 = oumega-H00;
t0 =(oumega_00)\H01';
t0_bar =(oumega_00)\H01;
T = t0;
T_bar = t0_bar;
ti = t0;
ti_bar = t0_bar;
ti_times = t0;
ti_bar_times = t0_bar;
I = eye(length(H00));
for i = 1:mu_max
    M1 = ti*ti_bar;
    M2 = ti_bar*ti;
    M3 = M1 + M2;
    ti = (I-M3)\ti^2;
    ti_bar = (I-M3)\ti_bar^2;
    T =T + ti_bar_times*ti;
    T_bar =T_bar + ti_times*ti_bar;
    %difference = ti +ti_bar;
    if norm(M3,'fro') < infinity_small
        %disp('reach_accuracy')
        break;
    end

    ti_times = ti_times*ti;
    ti_bar_times = ti_bar_times*ti_bar;

end

Gwl = inv(oumega_00-H01*T);
Gwb = inv((oumega_00- H01*T-H01'*T_bar));
Gwr = inv((oumega_00-H01'*T_bar));
end