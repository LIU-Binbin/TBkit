% function Arot = rotateinner(A,abc,RightorLeft,immproper,conjugate,antisymmetry)
% Arot =  rotateinner@Y_l__m(A,abc,RightorLeft,immproper,conjugate,antisymmetry);
% end

function Ak = rotateinner(A,abc,RightorLeft,immproper,conjugate,antisymmetry)
% rotateinner - Rotates the spherical harmonics using Wigner D-matrix
% Arguments:
% obj - Current object
% abc - Rotation angles [alpha, beta, gamma]
% RightorLeft - Flag for right or left rotation
% immproper, conjugate, antisymmetry - Optional parameters for symmetry handling

arguments
    A
    abc
    RightorLeft
    immproper = true;
    conjugate = false;
    antisymmetry = false;
end
alpha = (abc(1));
beta = (abc(2));
gamma = (abc(3));
A_L = Y_l__m(A.l);
Ak = A.HollowMe;
for i =  1:length(A_L)
    Ai = Y_lm(A_L(i));
    if A.l == 0
        Ai.coe = 1;
    else
        m1 = A.m;
        m2 = Ai.m;
        WignerD_single_element = (Y_l__m.d(A.l,m1,m2,beta));
        Ai.coe = Ai.coe*exp(1i*RightorLeft*m1*alpha)*WignerD_single_element*exp(1i*RightorLeft*m2*gamma);
        % for Y_l__m
        Ai.coe = conj(Ai.coe);
    end
    if immproper
        Ai.coe = (-1)^(A.l)*Ai.coe;
    end
    if ~conjugate
        Ai.coe = Ai.coe * A.coe;
    else
        Ai.coe = conj(Ai.coe) * A.coe; % check why the later cant use conj!
    end
    if Ai.coe ~= zeros(1,1,class(Ai.coe)) && ~isnan(Ai.coe)
        Ak = [Ak,Ai];
    end
end
end

