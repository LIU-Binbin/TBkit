function Ak = rotateinner(A, abc, RightorLeft, immproper, conjugate, antisymmetry)
arguments
A
abc
RightorLeft
immproper = false;
conjugate = false;
antisymmetry = false;
end
alpha = abc(1);
beta = abc(2);
gamma = abc(3);
A_L = Spin(A.J);
Ak = A.HollowMe;
for i = 1:length(A_L)
Ai = A_L(i);
if A.J == 0
Ai.coe = 1;
else
m1 = A.Jz;
m2 = Ai.Jz;
WignerD_single_element = Y_l__m.d(A.J, m1, m2, beta);
Ai.coe = Ai.coe * exp(1i * RightorLeft * m1 * alpha) * WignerD_single_element * exp(1i * RightorLeft * m2 * gamma);
end
if immproper
if A.parity == 1
Ai.coe = Ai.coe;
else
Ai.coe = -Ai.coe;
end
end
if ~conjugate
Ai.coe = Ai.coe * A.coe;
else
Ai = TimeRerversal(Ai);
Ai.coe = conj(Ai.coe * A.coe);
end
Ak = [Ak, Ai];
end
end
