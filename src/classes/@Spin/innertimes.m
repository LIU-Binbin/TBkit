function C = innertimes(A, B)
[Ergebnis, S, Sz] = CGM(A.J, B.J, A.Jz, B.Jz);
C = Spin(S, Sz, Ergebnis);
end
