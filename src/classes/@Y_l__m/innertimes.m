function C = innertimes(A,B)
if A.n ~= B.n
C = 0;
return;
end
l1 = A.l;l2 = B.l;
m1 = A.m;m2 = B.m;
M = m1+m2;
count = 0;
if A.symbolic ||B.symbolic
CoePublic = sqrt((2*l1+1)*(2*l2+1)/sym(pi))/2 ;
for L = abs(l1-l2):abs(l1+l2)
C1 = sym((-1)^L*sqrt(2*L+1));
W1 = Y_l__m.Wigner3j_sym([l1 l2 L],[m1 m2 -M]);
W2 = Y_l__m.Wigner3j_sym([l1 l2 L],[0 0 0]);
if W1~=sym(0) && W2~=sym(0)
count = count +1;
lL(count) = L;
COE(count) = C1*W1*W2;
end
end
else
CoePublic = sqrt((2*l1+1)*(2*l2+1)/pi)/2 ;
for L = abs(l1-l2):abs(l1+l2)
C1 = (-1)^L*sqrt(2*L+1);
W1 = Y_l__m.Wigner3j([l1 l2 L],[m1 m2 -M]);
W2 = Y_l__m.Wigner3j([l1 l2 L],[0 0 0]);
if W1~=0 && W2~=0
count = count +1;
lL(count) = L;
COE(count) = C1*W1*W2;
end
end
end
C  = Y_l__m(lL,M,CoePublic*COE,A.n);
end
